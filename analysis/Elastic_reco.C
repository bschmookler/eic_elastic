#include "Boost.h"

//------------------
void Elastic_reco(){

    //Define Style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(1.0);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);

    //Define histograms
    TH2 *h1_elec = new TH2D("h1_elec","Scattered electron true momentum vs. polar angle",100,110,170,100,4.6,6.8);
    h1_elec->GetXaxis()->SetTitle("#theta [deg]"); h1_elec->GetXaxis()->CenterTitle();
    h1_elec->GetYaxis()->SetTitle("p [GeV/c]"); h1_elec->GetYaxis()->CenterTitle();
   
    TH2 *h1_proton = new TH2D("h1_proton","Outgoing proton true momentum vs. polar angle",100,0,15,100,20,60);
    h1_proton->GetXaxis()->SetTitle("#theta [deg]"); h1_proton->GetXaxis()->CenterTitle();
    h1_proton->GetYaxis()->SetTitle("p [GeV/c]"); h1_proton->GetYaxis()->CenterTitle();
 
    TH2 *h2_elec = new TH2D("h2_elec","Scattered electron reconstructed momentum vs. polar angle",100,110,170,100,4.6,6.8);
    h2_elec->GetXaxis()->SetTitle("#theta [deg]"); h2_elec->GetXaxis()->CenterTitle();
    h2_elec->GetYaxis()->SetTitle("p [GeV/c]"); h2_elec->GetYaxis()->CenterTitle();

    TH2 *h2_proton = new TH2D("h2_proton","Outgoing proton reconstructed momentum vs. polar angle",100,0,15,100,20,60);
    h2_proton->GetXaxis()->SetTitle("#theta [deg]"); h2_proton->GetXaxis()->CenterTitle();
    h2_proton->GetYaxis()->SetTitle("p [GeV/c]"); h2_proton->GetYaxis()->CenterTitle();

    TH1 *h3a = new TH1D("h3a","Reconstructed proton and electron #phi balance: Lab frame",100,-600,600);
    h3a->SetLineWidth(2);h3a->SetLineColor(kBlue);
    h3a->GetXaxis()->SetTitle("|#phi_{prot.} - #phi_{elec.}| - #pi [mRad]");h3a->GetXaxis()->CenterTitle();

    TH1 *h3b = new TH1D("h3b","Reconstructed proton and electron P_{T} balance: Lab frame",100,0,2);
    h3b->SetLineWidth(2);h3b->SetLineColor(kBlue);
    h3b->GetXaxis()->SetTitle("P_{T,prot.} / P_{T,elec.}");h3b->GetXaxis()->CenterTitle();

    TH1 *h3c = new TH1D("h3c","Reconstructed total E - p_{z}: Lab frame",100,5,15);
    h3c->SetLineWidth(2);h3c->SetLineColor(kBlue);
    h3c->GetXaxis()->SetTitle("Total E - p_{z} [GeV]");h3c->GetXaxis()->CenterTitle();

    TH1 *h4a = new TH1D("h4a","Reconstructed proton and electron #phi balance: Colinear frame",100,-600,600);
    h4a->SetLineWidth(2);h4a->SetLineColor(kBlue);
    h4a->GetXaxis()->SetTitle("|#phi_{prot.} - #phi_{elec.}| - #pi [mRad]");h4a->GetXaxis()->CenterTitle();

    TH1 *h4b = new TH1D("h4b","Reconstructed proton and electron P_{T} balance: Colinear frame",100,0,2);
    h4b->SetLineWidth(2);h4b->SetLineColor(kBlue);
    h4b->GetXaxis()->SetTitle("P_{T,prot.} / P_{T,elec.}");h4b->GetXaxis()->CenterTitle();

    TH1 *h4c = new TH1D("h4c","Reconstructed total E - p_{z}: Colinear frame",100,5,15);
    h4c->SetLineWidth(2);h4c->SetLineColor(kBlue);
    h4c->GetXaxis()->SetTitle("Total E - p_{z} [GeV]");h4c->GetXaxis()->CenterTitle();

    TH2* h5a = new TH2D("h5a","Q^{2} reconstruction: electron method",100,0,40,100,0,40);
    h5a->GetXaxis()->SetTitle("Q^{2}_{true}");h5a->GetXaxis()->CenterTitle();
    h5a->GetYaxis()->SetTitle("Q^{2}_{elec.}");h5a->GetYaxis()->CenterTitle();

    TH2* h5b = new TH2D("h5b","Q^{2} reconstruction: J.B. method",100,0,40,100,0,40);
    h5b->GetXaxis()->SetTitle("Q^{2}_{true}");h5b->GetXaxis()->CenterTitle();
    h5b->GetYaxis()->SetTitle("Q^{2}_{J.B.}");h5b->GetYaxis()->CenterTitle();

    TH2* h5c = new TH2D("h5c","Q^{2} reconstruction: D.A. method",100,0,40,100,0,40);
    h5c->GetXaxis()->SetTitle("Q^{2}_{true}");h5c->GetXaxis()->CenterTitle();
    h5c->GetYaxis()->SetTitle("Q^{2}_{D.A.}");h5c->GetYaxis()->CenterTitle();

    TH2* h6a = new TH2D("h6a","x reconstruction: electron method",100,0,40,100,0,2);
    h6a->GetXaxis()->SetTitle("Q^{2}_{true}");h6a->GetXaxis()->CenterTitle();
    h6a->GetYaxis()->SetTitle("x_{elec.}");h6a->GetYaxis()->CenterTitle();

    TH2* h6b = new TH2D("h6b","x reconstruction: J.B. method",100,0,40,100,0,2);
    h6b->GetXaxis()->SetTitle("Q^{2}_{true}");h6b->GetXaxis()->CenterTitle();
    h6b->GetYaxis()->SetTitle("x_{J.B.}");h6b->GetYaxis()->CenterTitle();

    TH2* h6c = new TH2D("h6c","x reconstruction: D.A. method",100,0,40,100,0,2);
    h6c->GetXaxis()->SetTitle("Q^{2}_{true}");h6c->GetXaxis()->CenterTitle();
    h6c->GetYaxis()->SetTitle("x_{D.A.}");h6c->GetYaxis()->CenterTitle();

    //Works in EIC-shell
    //May give errors when run outside eic-shell
    TChain *tree = new TChain("events");
    tree->Add("input/eicrecon_*.root");

    //Print number of events
    int nevents = tree->GetEntries();
    cout<<"--------------"<<endl;
    cout<<"The total number of events to analyze is "<<nevents<<endl;
    cout<<"--------------"<<endl;

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<float> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<float> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<float> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass");
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");
  
    TTreeReaderArray<int> rec_pid(tr, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> rec_px(tr, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_py(tr, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_pz(tr, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_mass(tr, "ReconstructedChargedParticles.mass");

    //Other variables
    TLorentzVector gen_vec; //generated particle in lab frame
    TVector3 gen_vertex; //generated vertex in lab frame
    TLorentzVector rec_vec; //reconstructed particle in lab frame
    TLorentzVector rec_vec_co; //reconstructed particle boosted to colinear frame
    
    TLorentzVector e_beam_true; //Electron beam including event-by-event fluctuations
    TLorentzVector p_beam_true; //Proton beam including event-by-event fluctuations

    const double crossingAngle = -0.025;
    TLorentzVector e_beam(0.,0.,-5.,5.); //Electron beam for reconstruction (no event-by-event fluctuations)
    TLorentzVector p_beam(41.*sin(crossingAngle),0,41.*cos(crossingAngle),41.); //Proton beam for reconstruction (no event-by-event fluctuations)
    const double s_cm = 4.*5.*41.; //Square of center-of-mass energy
    
    TLorentzVector elec_gen; //Generated scattered elctron in lab frame
    TLorentzVector elec_rec; //Reconstructed scattered electron from tracking detector in lab frame

    bool found_elec; //Reconstruction status of scattered electron
    bool found_prot; //Reconstruction status of scattered proton
    TLorentzVector prot_rec; //Reconstructed proton from tracking detector in lab frame
    TLorentzVector elec_rec_co; //Reconstructed scattered electron boosted to colinear frame
    TLorentzVector prot_rec_co; //Reconstructed proton boosted to colinear frame
    
    int counter(0);

    //Loop over events
    while (tr.Next()) {
	
	if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	counter++;

        //Reset variables
        found_elec = false; found_prot = false;

        //Loop over generated particles, select primary electron, pi-plus, neutron
        for(int igen=0;igen<gen_status.GetSize();igen++){

                //Get Electron true beam momentum for event
                if(gen_status[igen]==4 && gen_pid[igen]==11){
			e_beam_true.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                } 
                //Get Proton true beam momentum for event
                if(gen_status[igen]==4 && gen_pid[igen]==2212){
			p_beam_true.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                } 
                //Get true scattered electron and outgoing proton
        	if(gen_status[igen]==1){
                	gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                	gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
			
			if(gen_pid[igen]==11){
                                h1_elec->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
                                elec_gen = gen_vec;
                        }
			if(gen_pid[igen]==2212) h1_proton->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
                }

        } //End loop over generated particles

        //Calculate true Q2 -- can use InclusiveKinematicsTruth collection instead
        //Using 4-vector algebra, so calculation will be same in lab and colinear frame
        auto Q2_true = -1. * ((e_beam_true - elec_gen).Mag2());

	//Loop over reconstructed tracks for electron and pion
	for(int irec=0;irec<rec_pid.GetSize();irec++){
		rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);
                rec_vec_co = apply_boost(e_beam,p_beam,rec_vec);

		if(rec_pid[irec]==11){
			h2_elec->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
                        elec_rec = rec_vec;
                        elec_rec_co = rec_vec_co;
                        found_elec = true;
		}

                if(rec_pid[irec]==2212){
			h2_proton->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
                        prot_rec = rec_vec;
                        prot_rec_co = rec_vec_co;
                        found_prot = true;
		}

	}//End loop over reconstructed particles

        //Calculate reconstructed Q2 and x values
        //---
        //Electron method -- use 4-vector algebra, so can do calculation either in lab or colinear frame
        if(found_elec){
                auto q_elec = e_beam - elec_rec;
                auto Q2_elec = -1.*q_elec*q_elec;
                auto x_elec = Q2_elec / (2.*p_beam*q_elec);

                h5a->Fill(Q2_true,Q2_elec);
                h6a->Fill(Q2_true,x_elec);
        }

        //J.B. method -- reconstruction using reconstructed proton (hadronic final state)
        //We do the reconstruction in the colinear frame. 
        //Note the electron and proton beam energies are (almost) identical in both frames
        if(found_prot){
                auto sigma_h = prot_rec_co.E() - prot_rec_co.Pz();
                auto pt_h = prot_rec_co.Pt();

                auto y_jb = sigma_h / (2.*e_beam.E());
                auto Q2_jb = (pt_h*pt_h) / (1. - y_jb);
                auto x_jb = Q2_jb / (s_cm*y_jb);

                h5b->Fill(Q2_true,Q2_jb);
                h6b->Fill(Q2_true,x_jb);
        }

        //D.A method -- reconstruction using electron and proton polar angles
        //We do the reconstruction in the colinear frame. 
        //Note the electron and proton beam energies are (almost) identical in both frames
        if( found_elec && found_prot ){
                auto tan_tpo2 = std::tan( prot_rec_co.Theta()/2. );
                auto tan_teo2 = std::tan( elec_rec_co.Theta()/2. );
                auto cot_teo2 = 1./tan_teo2;

                auto y_da = (tan_tpo2) / (tan_tpo2 + tan_teo2);
                auto Q2_da = 4. * std::pow(e_beam.E(),2) * (cot_teo2) / (tan_tpo2 + tan_teo2);
                auto x_da = Q2_da / (s_cm*y_da);

                h5c->Fill(Q2_true,Q2_da);
                h6c->Fill(Q2_true,x_da);
        }
        //---

        //Fill additional kinematic histograms
        //Require reconstruction of both the scattered electron and the outgoing proton
        //---
        if(found_elec && found_prot){
                //Lab frame
                auto phi_diff = std::abs(prot_rec.Phi() - elec_rec.Phi()) - M_PI;
                auto Empz_tot = prot_rec.E() - prot_rec.Pz() + elec_rec.E() - elec_rec.Pz();

                h3a->Fill(1000.*phi_diff);
                h3b->Fill(prot_rec.Pt() / elec_rec.Pt());
                h3c->Fill(Empz_tot);

                //Colinear frame
                auto phi_diff_co = std::abs(prot_rec_co.Phi() - elec_rec_co.Phi()) - M_PI;
                auto Empz_tot_co = prot_rec_co.E() - prot_rec_co.Pz() + elec_rec_co.E() - elec_rec_co.Pz();

                h4a->Fill(1000.*phi_diff_co);
                h4b->Fill(prot_rec_co.Pt() / elec_rec_co.Pt());
                h4c->Fill(Empz_tot_co);
        }
        //---

     } //End loop over events

    //Make plots
    TCanvas *c1a = new TCanvas("c1a");
    h1_elec->Draw("colz");

    TLatex *tex1 = new TLatex();
    tex1->SetTextSize(0.035);
    tex1->DrawLatex(140,6.65,"Elastic generator: e + p #rightarrow e + p");
    tex1->DrawLatex(140,6.55,"5x41 GeV: Q^{2}_{gen.} > 5 GeV^{2}/c^{2}");
    tex1->DrawLatex(140,6.45,"Beam-effects afterburner applied");
    tex1->DrawLatex(140,6.25,"100k events generated: #intL = 7.3 fb^{-1}");
    TCanvas *c1b = new TCanvas("c1b");
    h1_proton->Draw("colz");

    TCanvas *c2a = new TCanvas("c2a");
    h2_elec->Draw("colz");

    TCanvas *c2b = new TCanvas("c2b");
    h2_proton->Draw("colz");

    TCanvas *c3a = new TCanvas("c3a");
    h3a->Draw("");
    
    TCanvas *c3b = new TCanvas("c3b");
    h3b->Draw("");

    TCanvas *c3c = new TCanvas("c3c");
    h3c->Draw("");

    TCanvas *c4a = new TCanvas("c4a");
    h4a->Draw("");
    
    TCanvas *c4b = new TCanvas("c4b");
    h4b->Draw("");

    TCanvas *c4c = new TCanvas("c4c");
    h4c->Draw("");

    TCanvas *c5a = new TCanvas("c5a");
    c5a->SetLogz();
    h5a->Draw("colz");

    TCanvas *c5b = new TCanvas("c5b");
    c5b->SetLogz();
    h5b->Draw("colz");

    TCanvas *c5c = new TCanvas("c5c");
    c5c->SetLogz();
    h5c->Draw("colz");

    TCanvas *c6a = new TCanvas("c6a");
    c6a->SetLogz();
    h6a->Draw("colz");

    TCanvas *c6b = new TCanvas("c6b");
    c6b->SetLogz();
    h6b->Draw("colz");

    TCanvas *c6c = new TCanvas("c6c");
    c6c->SetLogz();
    h6c->Draw("colz");

    //Print plots to file
    c1a->Print("plots/Elastic_reco.pdf[");
    c1a->Print("plots/Elastic_reco.pdf");
    c1b->Print("plots/Elastic_reco.pdf");
    c2a->Print("plots/Elastic_reco.pdf");
    c2b->Print("plots/Elastic_reco.pdf");
    c3a->Print("plots/Elastic_reco.pdf");
    c3b->Print("plots/Elastic_reco.pdf");
    c3c->Print("plots/Elastic_reco.pdf");
    c4a->Print("plots/Elastic_reco.pdf");
    c4b->Print("plots/Elastic_reco.pdf");
    c4c->Print("plots/Elastic_reco.pdf");
    c5a->Print("plots/Elastic_reco.pdf");
    c5b->Print("plots/Elastic_reco.pdf");
    c5c->Print("plots/Elastic_reco.pdf");
    c6a->Print("plots/Elastic_reco.pdf");
    c6b->Print("plots/Elastic_reco.pdf");
    c6c->Print("plots/Elastic_reco.pdf");
    c6c->Print("plots/Elastic_reco.pdf]");

}
