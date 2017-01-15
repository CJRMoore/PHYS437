std::vector<double> angles;
std::vector<double> prob;
TGraph *g;
TH1D *h1;
fitgaussian(std::string file, std::string title="", std::string xtitle="", std::string ytitle=""){
    angles.resize(0);
    prob.resize(0);
    std::ifstream f(file.c_str());;


    while (f.good()){
        double a, p;
        f >> a >> p;
        angles.push_back(a);
        prob.push_back(p);
    }

    if (angles.back()<angles[0]) h1 = new TH1D("h1","h1",angles.size(),angles.back(),angles[0]);
    else h1 = new TH1D("h1","h1",angles.size(),angles[0],angles.back());
    for (int i=0; i<angles.size(); i++){
        h1->SetBinContent(h1->FindBin(angles[i]),prob[i]);
    }
    h1->Draw();
    if (title!="") h1->SetTitle(title.c_str());
    if (xtitle!="") h1->SetXTitle(xtitle.c_str());
    if (ytitle!="") h1->SetYTitle(ytitle.c_str());
    h1->Draw();
    //g = new TGraph(angles.size(),&angles[0],&prob[0]);
    //g->Draw();
}
