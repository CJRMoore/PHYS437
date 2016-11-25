#include <iostream>
#include "Particles/Molecule.h"
#include "Event/Field.h"
#include "Event/Event.h"
#include "Event/Event.cxx"

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"

int main(int argc, char** argv){
    int nIterations = 1;
    std::string outFile = "output.root";
    std::string openMethod = "RECREATE";
    if (argc==3) {
        nIterations = atoi(argv[1]);
        outFile = argv[2];
    }
    else if (argc==4) {
        nIterations = atoi(argv[1]);
        outFile = argv[2];
        openMethod = argv[3];
    }
    else if (argc==2) nIterations = atoi(argv[1]);

    std::vector<double> bondlengths(2,0);
    double bondangle = 0;
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
        tree->Branch("bondL",&bondlengths);
        tree->Branch("bondA",&bondangle);
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

        std::vector<double> *blp = &bondlengths;
        double *bap = &bondangle;
        tree->SetBranchAddress("mass",&mp);
        tree->SetBranchAddress("charge",&cp);
        tree->SetBranchAddress("x",&xp);
        tree->SetBranchAddress("y",&yp);
        tree->SetBranchAddress("tof",&tofp);
        tree->SetBranchAddress("px",&pxp);
        tree->SetBranchAddress("py",&pyp);
        tree->SetBranchAddress("pz",&pzp);
        tree->SetBranchAddress("bondL",&blp);
        tree->SetBranchAddress("bondA",&bap);
    }

    Field *f = new Field();

    for (int i=0; i<nIterations; i++){
        //std::cout << "\n\nNEW RUN" << std::endl;
//        std::cout << "================Initializing!================\n";
        std::cerr << "\rProgress: " << i+1 << "\tInitializing Molecules" << std::flush;
        Molecule *m = new Molecule();
        bondangle = m->GetAngle();
        bondlengths = m->GetBondLengths();
        m->Rotate();
        m->Ionize();

        /*for (int j=0; j<m->GetNatoms(); j++){
            Atom *atom = m->GetAtom(j);
            Eigen::Vector3d translate = atom->GetPosition();
            translate[2] += Zinitial;
            atom->SetPosition(translate);
        } */

        EventHandler *e = new EventHandler(f, m);
        std::cerr << "\rProgress: " << i+1 << "\tExplosion             " << std::flush;
        double ExplosionTime = e->Run(0);

        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            mass[j] = atom->GetMass();
            charge[j] = atom->GetTotalCharge();
            Eigen::Vector3d vel = atom->GetVelocity();
            px[j] = vel[0] * mass[j];
            py[j] = vel[1] * mass[j];
            pz[j] = vel[2] * mass[j];

            //Eigen::Vector3d translate = atom->GetPosition();
            //translate[2] += Zinitial;
            //atom->SetPosition(translate);
        }


//        std::cout << "Finished explosion in " << ExplosionTime << " seconds.\n";
        std::cerr << "\rProgress: " << i+1 << "\tExtraction            " << std::flush;
        double ToF = e->Run(1);



        bool failed=false;
        for (int j=0; j<m->GetNatoms(); j++){
            Atom* atom = m->GetAtom(j);
            x[j] = atom->GetPosition()[0];
            y[j] = atom->GetPosition()[1];
            tof[j] = atom->GetTimeOfFlight();
            if (fabs(x[j]+1)<=1e-12 || fabs(y[j]+1)<=1e-12) failed=true;
        }
        
        if (!failed){
            tree->Fill();
            //file->Write();
        }
        //else i--;
        delete m;
        delete e;
    }
    std::cout << std::endl;
    file->Write();
    //}
    delete f;

    delete tree;
    delete file;
}
