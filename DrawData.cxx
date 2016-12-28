TCanvas *c1 = 0;
TPad *p1 = 0;

DrawData(std::string file, int particle = -1, std::string buffer = ""){
    TFile *f = new TFile(file.c_str());
    if (!f) {
        std::cerr << "Missing TTree\n";
        return;
    }
      
    TTree *t = f->Get("data");
    if (!t) {
        std::cerr << "Missing TTree\n";
        return;
    }

    if (c1==0) {
        c1 = new TCanvas("c1","c1",100,100,800,600);
        c1->Draw();
    }

    //Row 1
    if (p1==0){
        TPad *p1 = new TPad("p1","p1",0.01,0.51,0.32,0.99);
        TPad *p2 = new TPad("p2","p2",0.33,0.51,0.65,0.99);
        TPad *p3 = new TPad("p3","p3",0.66,0.51,0.99,0.99);
    //Row 2
        TPad *p4 = new TPad("p4","p4",0.01,0.01,0.32,0.49);
        TPad *p5 = new TPad("p5","p5",0.33,0.01,0.65,0.49);
        TPad *p6 = new TPad("p6","p6",0.66,0.01,0.99,0.49);

        p1->Draw();
        p2->Draw();
        p3->Draw();
        p4->Draw();
        p5->Draw();
        p6->Draw();
    }

    // x position
    p1->cd();
    p1->SetTitle("X position (m)");
    char buf[100];
    sprintf(buf,"x");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());

    // y position
    p2->cd();
    p2->SetTitle("Y position (m)");
    sprintf(buf,"y");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());

    // tof
    p3->cd();
    p3->SetTitle("Time of Flight (s)");
    sprintf(buf,"tof");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());

    // px
    p4->cd();
    p4->SetTitle("X momentum");
    sprintf(buf,"px");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());

    // py
    p5->cd();
    p5->SetTitle("Y momentum");
    sprintf(buf,"py");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());

    // pz
    p6->cd();
    p6->SetTitle("Z momentum");
    sprintf(buf,"pz");
    if (particle!=-1) sprintf(buf,"%s[%i]",buf,particle);
    t->Draw(buf,buffer.c_str());
}
