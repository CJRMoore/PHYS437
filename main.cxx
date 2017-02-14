#include <iostream>
#include <omp.h>
#include <fstream>
#include <sstream>

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

    std::cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << " %\r";
    std::cerr.flush();
}


int main(int argc, char** argv){
    int nIterations = 1;
    std::string outFile = "";
    std::string openMethod = "RECREATE";
    std::string fieldFile = "geo_default.txt";
    int corecount = omp_get_num_threads();
    std::vector<int> ionizations;
    std::string seedFile = "";
    std::string atomSpec = "";
    bool varyField = false;

    for (int i=1; i<argc; i++){
        if (strcmp(argv[i],"-n")==0)        nIterations = atoi(argv[i+1]);
        else if (strcmp(argv[i],"-c")==0)     corecount = atoi(argv[i+1]);
        else if (strcmp(argv[i],"-o")==0)       outFile = argv[i+1];
        else if (strcmp(argv[i],"-m")==0)    openMethod = argv[i+1];
        else if (strcmp(argv[i],"-f")==0)     fieldFile = argv[i+1];
        else if (strcmp(argv[i],"-s")==0)      seedFile = argv[i+1];
        else if (strcmp(argv[i],"-a")==0)      atomSpec = argv[i+1];
        else if (strcmp(argv[i],"-v")==0)     varyField = true;
        else if (strcmp(argv[i],"-I")==0){
            int inputFromCL=0;
            int CLnum=1;
            std::istringstream ss;
            while (true){
                ss.str(argv[i+CLnum]);
                ss >> inputFromCL;
                if (ss.fail()) break;
                ionizations.push_back(inputFromCL);
                CLnum++;
            }
            /*I1 = atoi(argv[i+1]);
            I2 = atoi(argv[i+2]);
            I3 = atoi(argv[i+3]);
            i += 2;*/
            i += CLnum-1;
        }
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
    
    std::vector<unsigned int> seeds(0);
    if (seedFile!=""){
        if (nIterations!=1) std::cerr << "Overriding number of iterations (given seed file)\n";
        std::ifstream sfile(seedFile.c_str());
        std::string line;
        unsigned int seed;

        while (std::getline(sfile,line)){
            std::stringstream ss(line);
            ss >> seed;
            seeds.push_back(seed);
        }
        nIterations = seeds.size();
    }

    std::vector<double> mass(3,0);
    std::vector<double> charge(3,0);
    std::vector<double> x(3,0);
    std::vector<double> y(3,0);
    std::vector<double> x0(3,0);
    std::vector<double> y0(3,0);
    std::vector<double> z0(3,0);
    std::vector<double> tof(3,0);
    std::vector<double> px(3,0);
    std::vector<double> py(3,0);
    std::vector<double> pz(3,0);
    double kinetic=0;
    unsigned int Bseed = 0;
    double samplefield = 0;

    gROOT->ProcessLine("#include <vector>");

    TFile *file = new TFile(outFile.c_str(),openMethod.c_str());
    TTree *tree = (TTree*)file->Get("data");
    if (tree==0) {
        tree = new TTree("data","OCS explosion data");
        tree->Branch("mass",&mass);
        tree->Branch("charge",&charge);
        tree->Branch("x",&x);
        tree->Branch("y",&y);
        tree->Branch("x0",&x0);
        tree->Branch("y0",&y0);
        tree->Branch("z0",&z0);
        tree->Branch("tof",&tof);
        tree->Branch("px",&px);
        tree->Branch("py",&py);
        tree->Branch("pz",&pz);
        tree->Branch("ke",&kinetic);
        tree->Branch("seed",&Bseed);
        if (varyField) tree->Branch("field",&samplefield);
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
        // Make thread safe
            std::vector<double> tmass(3,0);
            std::vector<double> tcharge(3,0);
            std::vector<double> tx(3,0);
            std::vector<double> ty(3,0);
            std::vector<double> tx0(3,0);
            std::vector<double> ty0(3,0);
            std::vector<double> tz0(3,0);
            std::vector<double> ttof(3,0);
            std::vector<double> tpx(3,0);
            std::vector<double> tpy(3,0);
            std::vector<double> tpz(3,0);
            double tkinetic=0;
            unsigned int tBseed = 0;
            double tsamplefield = 0;

        #pragma omp critical
        {
            prognow++;
            ProgressBar((double) prognow, (double) nIterations);
        }
        
        Molecule *m;
        #pragma omp critical
        {
            if (seeds.size()>0) {
                tBseed = seeds[i];
                m = new Molecule(atomSpec,seeds[i]);
                m->Rotate(seeds[i]);
                m->GenerateVelocity(seeds[i]);
            }
            else {
                m = new Molecule(atomSpec);
                m->Rotate();
                m->GenerateVelocity();
            }
        }
        m->Ionize(ionizations);

        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            Eigen::Vector3d pos = atom->GetPosition();
            tx0[j] = pos[0];
            ty0[j] = pos[1];
            tz0[j] = pos[2];
        }
        
        EventHandler *e = new EventHandler(f, m);
        if (atomSpec!="D2") e->Run(0);

        tkinetic = m->GetKE();

        double Zinitial = 0.5 * (0.5 * (89.61+92.91) + 0.5 * (82.51+85.81)) * 1e-3;
        // for specific D2 test requested by Benji
        /*if (atomSpec=="D2") {
            Zinitial -= i * (Zinitial/(2*nIterations));
            tz0[0] = Zinitial;
        }*/
        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            tmass[j] = atom->GetMass();
            tcharge[j] = atom->GetTotalCharge();
            Eigen::Vector3d vel = atom->GetVelocity();
            tpx[j] = vel[0] * tmass[j];
            tpy[j] = vel[1] * tmass[j];
            tpz[j] = vel[2] * tmass[j];

            Eigen::Vector3d pos = atom->GetPosition();
            pos[2] += Zinitial;
            atom->SetPosition(pos);
        }

        if (varyField) e->Run(i+1);
        else e->Run(1);

        bool failed=false;
        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            tx[j] = atom->GetPosition()[0];
            ty[j] = atom->GetPosition()[1];
            ttof[j] = atom->GetTimeOfFlight();
            if (fabs(tx[j]+1)<=1e-12 || fabs(ty[j]+1)<=1e-12) failed=true;
        }

        #pragma omp critical
        {
            if (!failed) {
                for (int j=0; j<m->GetNatoms(); j++){
                    mass[j] = tmass[j];
                    charge[j] = tcharge[j];
                    x[j] = tx[j];
                    y[j] = ty[j];
                    x0[j] = tx0[j];
                    y0[j] = ty0[j];
                    z0[j] = tz0[j];
                    tof[j] = ttof[j];
                    px[j] = tpx[j];
                    py[j] = tpy[j];
                    pz[j] = tpz[j];
                }
                kinetic = tkinetic;
                Bseed = tBseed;
                if (varyField){
                    // Should only be used with static field
                    Eigen::Vector3d tempV(0,0,1e-3);
                    tempV = f->GetFieldAtPosition(tempV);
//                    std::cout << tempV[2] << std::endl;
                    samplefield = tempV[2] * (1001.-i)/1000;
                }
                tree->Fill();
            }
        }
        delete m;
        delete e;
    }

    std::cerr << std::endl;
    file->Write();
    delete f;

    delete tree;
    delete file;
}
