#include <iostream>
#include <omp.h>
#include "Particles/Molecule.h"
#include "Event/Field.h"
#include "Event/Event.h"
#include "Event/Event.cxx"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"


void ProgressBar(double now, double total){
    float progress = now/total;
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}


int main(int argc, char** argv){
    int nIterations = 1;
    std::string outFile = "";
    std::string openMethod = "RECREATE";
    std::string fieldFile = "geo_default.txt";
    int corecount = omp_get_num_threads();

    for (int i=1; i<argc; i++){
        if (strcmp(argv[i],"-n")==0)        nIterations = atoi(argv[i+1]);
        else if (strcmp(argv[i],"-c")==0)     corecount = atoi(argv[i+1]);
        else if (strcmp(argv[i],"-o")==0)       outFile = argv[i+1];
        else if (strcmp(argv[i],"-m")==0)    openMethod = argv[i+1];
        else if (strcmp(argv[i],"-f")==0)     fieldFile = argv[i+1];
        else continue;
        i++;
    }
    if (outFile=="") {
        std::cerr << "Missing output file; use argument -o <output file>\n";
        std::cerr << "Use default output file? (output.root)\t[Y/n]\n";
        char filebuf[100];
        std::cin >> filebuf;
        if (strcmp(filebuf,"n")==0 || strcmp(filebuf,"N")==0){
            std::cerr << "Breaking...\n";
            return 0;
        }
        else if (strcmp(filebuf,"y")==0 || strcmp(filebuf,"Y")==0){
            std::cerr << "Using default output file\n";
            outFile = "output.root";
        }
        else {
            std::cerr << "Unknown input: " << filebuf << std::endl;
            std::cerr << "Breaking...\n";
            return 0;
        }
    }

    omp_set_num_threads(corecount);
    

/*    if (argc==3) {
        nIterations = atoi(argv[1]);
        outFile = argv[2];
    }
    else if (argc==4) {
        nIterations = atoi(argv[1]);
        outFile = argv[2];
        openMethod = argv[3];
    }
    else if (argc==2) nIterations = atoi(argv[1]);*/

    std::vector<double> mass(3,0);
    std::vector<double> charge(3,0);
    std::vector<double> x(3,0);
    std::vector<double> y(3,0);
    std::vector<double> tof(3,0);
    std::vector<double> px(3,0);
    std::vector<double> py(3,0);
    std::vector<double> pz(3,0);

    gROOT->ProcessLine("#include <vector>");

    TFile *file = new TFile(outFile.c_str(),openMethod.c_str());
    TTree *tree = (TTree*)file->Get("data");
    if (tree==0) {
        tree = new TTree("data","OCS explosion data");
        tree->Branch("mass",&mass);
        tree->Branch("charge",&charge);
        tree->Branch("x",&x);
        tree->Branch("y",&y);
        tree->Branch("tof",&tof);
        tree->Branch("px",&px);
        tree->Branch("py",&py);
        tree->Branch("pz",&pz);
    }
    else{
        std::vector<double> *mp = &mass;
        std::vector<double> *cp = &charge;
        std::vector<double> *xp = &x;
        std::vector<double> *yp = &y;
        std::vector<double> *tofp = &tof;
        std::vector<double> *pxp = &px;
        std::vector<double> *pyp = &py;
        std::vector<double> *pzp = &pz;

        tree->SetBranchAddress("mass",&mp);
        tree->SetBranchAddress("charge",&cp);
        tree->SetBranchAddress("x",&xp);
        tree->SetBranchAddress("y",&yp);
        tree->SetBranchAddress("tof",&tofp);
        tree->SetBranchAddress("px",&pxp);
        tree->SetBranchAddress("py",&pyp);
        tree->SetBranchAddress("pz",&pzp);
    }

    Field *f = new Field(fieldFile);

    int prognow = 0;
    
    #pragma omp parallel for
    for (int i=0; i<nIterations; i++){
        #pragma omp critical
        {
            prognow++;
            ProgressBar((double) prognow, (double) nIterations);
        }
        
        Molecule *m;
        #pragma omp critical
        {
            m = new Molecule();
        }
        m->Rotate();
        m->Ionize();

        EventHandler *e = new EventHandler(f, m);
        double ExplosionTime = e->Run(0);

        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            mass[j] = atom->GetMass();
            charge[j] = atom->GetTotalCharge();
            Eigen::Vector3d vel = atom->GetVelocity();
            px[j] = vel[0] * mass[j];
            py[j] = vel[1] * mass[j];
            pz[j] = vel[2] * mass[j];
        }

        double ToF = e->Run(1);

        bool failed=false;
        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            x[j] = atom->GetPosition()[0];
            y[j] = atom->GetPosition()[1];
            tof[j] = atom->GetTimeOfFlight();
            if (fabs(x[j]+1)<=1e-12 || fabs(y[j]+1)<=1e-12) failed=true;
        }

        #pragma omp critical
        {
            if (!failed) tree->Fill();
        }
        delete m;
        delete e;
    }

    std::cout << std::endl;
    file->Write();
    delete f;

    delete tree;
    delete file;
}
