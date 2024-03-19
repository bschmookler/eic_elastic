#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

using namespace HepMC3;

//Globals
TLorentzVector Ei,Pi;
TLorentzVector Ei_rest,Pi_rest;
TLorentzVector Ef_rest,Pf_rest;
TLorentzVector Ef,Pf;
TLorentzVector q4vec;
const Double_t Mp(0.9383);
const Double_t Me(0.511E-3);
const Double_t alpha(1./137.036);
const Double_t mu_p(2.79285);

//Initial Generation Ranges
Double_t theta_low = 1E-6;
Double_t theta_hi = (TMath::Pi()/2.); 
Double_t phi_low = -TMath::Pi(); 
Double_t phi_hi = TMath::Pi();
Double_t cos_theta_low(1.);
Double_t cos_theta_hi(0.); 
Double_t Q2_low(0.);
Double_t Q2_hi(0.);

//-------------------------
//Dipole Form Factor GD
Double_t GD(Double_t Q2){
    return 1./(TMath::Power(1. + Q2/0.71,2));
}

//-------------------------
//GEp from Kelly Parameterization
Double_t GEp(Double_t Q2){

    Double_t gep_par[4] = {-.299, 11.11, 14.11, 15.7};
    Double_t a1 = gep_par[0];
    Double_t b1 = gep_par[1];
    Double_t b2 = gep_par[2];
    Double_t b3 = gep_par[3];

    Double_t tau = Q2/(4.*Mp*Mp);

    Double_t GE = (1. + a1*tau)/(1.+b1*tau+b2*pow(tau,2.)+b3*pow(tau,3.));
    return GE;
}

//-------------------------
//GMp from Kelly Parameterization
Double_t GMp(Double_t Q2){

    Double_t gmp_par[4] = {.081,11.15,18.45,5.31};
    Double_t a1 = gmp_par[0];
    Double_t b1 = gmp_par[1];
    Double_t b2 = gmp_par[2];
    Double_t b3 = gmp_par[3];

    Double_t tau = Q2/(4.*Mp*Mp);

    Double_t GM = mu_p*(1. + a1*tau)/(1.+b1*tau+b2*pow(tau,2.)+b3*pow(tau,3.));
    return GM;
}

//-------------------------
Double_t diff_cs(Int_t nDim, Double_t *Xarg){ //Will return cs in fb/dtheta*dphi in proton rest frame

    //Cross Section to return
    Double_t cs(0);
    Double_t cs_ns(0);
    Double_t ff_term(0);

    //Sampled variables
    Double_t theta = theta_low + Xarg[0]*(theta_hi - theta_low); //This theta will be relative to electron direction
    Double_t phi = phi_low + Xarg[1]*(phi_hi - phi_low);

    //Calculate Variables
    Double_t Eei = Ei_rest.E();
    Double_t Eef = ( Mp*Eei ) /( Mp + Eei*(1.-TMath::Cos(theta)) ); //Final Electron Energy
    Double_t Q2 = 4.*Eei*Eef*TMath::Sin(theta/2.)*TMath::Sin(theta/2.);
    Double_t tau = Q2/(4.*Mp*Mp);
    
    //Double_t GE = GD(Q2);  //Dipole FF
    //Double_t GM = mu_p*GE; //''
    Double_t GE = GEp(Q2);
    Double_t GM = GMp(Q2);

    cs_ns = (alpha*alpha)/(4.*Eei*Eei*TMath::Power(TMath::Sin(theta/2.),4))*Eef/Eei;
    ff_term = ( ((GE*GE)+(tau*GM*GM))*TMath::Power(TMath::Cos(theta/2.),2) / (1.+tau) );
    ff_term+= ( 2*tau*GM*GM*TMath::Power(TMath::Sin(theta/2.),2) );
    
    cs = cs_ns*ff_term*TMath::Sin(theta)*(0.3894)*(1E12);
    return cs;
}

//-------------------------
Double_t diff_cs_Q2(Double_t *x, Double_t *par){ //Will return cs in fb/GeV^2
  
  Double_t Q2 = x[0];
  Double_t tau = Q2/(4.*Mp*Mp);
  
  Double_t s_cm = par[0];

  Double_t GE = GEp(Q2);
  Double_t GM = GMp(Q2);

  Double_t y = Q2/( (par[0]) - (Mp*Mp) ); //Q2 = (s-Mp^2)*x*y; x=1 for elastic

  Double_t dsig_dQ2(0);
  
  if(y<=1){
    dsig_dQ2 = 4.*TMath::Pi()*(TMath::Power(alpha,2))/(TMath::Power(Q2,2));
    dsig_dQ2 = dsig_dQ2*( (((GE*GE + tau*GM*GM)/(1.+tau))*(1. - y -  (Mp*Mp*y*y/x[0]))) + 
      ((y*y/2.)*GM*GM) );
    dsig_dQ2 = dsig_dQ2*(0.3894)*(1E12);
  }
  else
    dsig_dQ2 = 0;

  return dsig_dQ2;

}

//-------------------------
void elas_gen_hepmc(){

    //Read In Beam Energies
    Int_t energy_set(0);

    cout<<"Enter Beam Energy Combination:"<<endl;
    cout<<"1) e:5 GeV; p:41 GeV"<<endl; //EIC Energy 1
    cout<<"2) e:5 GeV; p:100 GeV"<<endl; //EIC Energy 2
    cout<<"3) e:10 GeV; p:100 GeV"<<endl; //EIC Energy 3
    cout<<"4) e:10 GeV; p:275 GeV"<<endl; //EIC Energy 4
    cout<<"5) e:18 GeV; p:275 GeV"<<endl; //EIC Energy 5
    cin>>energy_set;

    TString kin;

    if(energy_set==1){
        cout<<"Doing 5x41 GeV combination"<<endl;
        Ei.SetXYZM(0,0,-5,Me);
        Pi.SetXYZM(0,0,41,Mp);

        kin = "5_41";
    }
    else if(energy_set==2){
        cout<<"Doing 5x100 GeV combination"<<endl;
        Ei.SetXYZM(0,0,-5,Me);
        Pi.SetXYZM(0,0,100,Mp);

        kin = "5_100";
    }
    else if(energy_set==3){
        cout<<"Doing 10x100 GeV combination"<<endl;
        Ei.SetXYZM(0,0,-10,Me);
        Pi.SetXYZM(0,0,100,Mp);

        kin = "10_100";
    }
    else if(energy_set==4){
        cout<<"Doing 10x275 GeV combination"<<endl;
        Ei.SetXYZM(0,0,-10,Me);
        Pi.SetXYZM(0,0,275,Mp);

        kin = "10_275";
    }
    else if(energy_set==5){
        cout<<"Doing 18x275 GeV combination"<<endl;
        Ei.SetXYZM(0,0,-18,Me);
        Pi.SetXYZM(0,0,275,Mp);

        kin = "18_275";
    }
  
    Double_t s_cm = (Ei+Pi).M2();

    //Boost
    TVector3 boost_vec = Pi.BoostVector();
    Ei_rest = Ei; Ei_rest.Boost(-boost_vec);
    Pi_rest = Pi; Pi_rest.Boost(-boost_vec);

    cout<<endl<<"Incoming Electron 4 Momentum in Lab Frame:"<<endl;
    Ei.Print();
    cout<<"Incoming Proton 4 Momentum in Lab Frame:"<<endl;
    Pi.Print();

    cout<<"Incoming Electron 4 Momentum in Proton Rest Frame:"<<endl;
    Ei_rest.Print();
    cout<<"Incoming Proton 4 Momentum in Proton Rest Frame:"<<endl;
    Pi_rest.Print();

    //Get Q2 Limits
    cout<<endl<<"Enter Q2 lower limit (in GeV^2):"<<endl;
    cin>>Q2_low;

    cout<<endl<<"Enter Q2 higher limit (in GeV^2):"<<endl;
    cin>>Q2_hi;

    //Calculate theta limits
    Double_t Ep_low = Ei_rest.E() - ( Q2_low/(2.*Mp) );
    cos_theta_low = 1. - ( Q2_low /(2*Ei_rest.E()*Ep_low) );
    theta_low = acos(cos_theta_low);

    Double_t Ep_hi = Ei_rest.E() - ( Q2_hi/(2.*Mp) );
    cos_theta_hi = 1. - ( Q2_hi /(2*Ei_rest.E()*Ep_hi) );
    theta_hi = acos(cos_theta_hi);

    //Calculate Integrated Cross Section
    TF1 *f1 = new TF1("f1",diff_cs_Q2,Q2_low,Q2_hi,1);
    f1->SetParameter(0,s_cm);
    Double_t CS_Int_Q2 = f1->Integral(Q2_low,Q2_hi); //In fb

    //Choose whether to use Foam
    Int_t do_foam(0);
    cout<<endl<<"Would you like to use the Foam of Uniform Generator:"<<endl;
    cout<<"0) Uniform"<<endl;
    cout<<"1) Foam"<<endl;
    cin>>do_foam;

    //Number of events to generate
    Int_t nevents;
    cout<<endl<<"How many events would you like to generate:"<<endl;
    cin>>nevents;

    Double_t Eei = Ei_rest.E();
    Double_t Eef(0);

    const Int_t num_dim = 2;
    Double_t *MCvect = new Double_t[num_dim];
    Double_t MCwt(0);
    Double_t cs_sum(0);

    Double_t theta_lab(0),theta_event(0),phi_event(0);
    Double_t Q2(0);

    //Open output HepMC3 file
    WriterAscii hepmc_output( Form("elas_gen_%s.hepmc",kin.Data()) );
    GenEvent evt(Units::GEV, Units::MM);

    //Generate Events
    TRandom3 *rand = new TRandom3(0);

    TFoam *Foam = new TFoam("Foam");
    if(do_foam==1){
        Foam->SetkDim(num_dim);
        Foam->SetRhoInt(diff_cs);
        Foam->SetPseRan(rand);
        Foam->Initialize();
    }

    //Loop
    for(Int_t iEvent=0;iEvent<nevents;iEvent++){
        if(iEvent%100000==0) cout<<"Generated Event "<<iEvent<<"!"<<endl;

        //Set the event number
        evt.set_event_number(iEvent);

        if(do_foam==1){
            Foam->MakeEvent();
            Foam->GetMCvect(MCvect);
            MCwt = Foam->GetMCwt(); //Gives 1
        
            theta_lab = theta_low + MCvect[0]*(theta_hi-theta_low);
            theta_event = TMath::Pi() - theta_lab; //Since initial electron is in -z direction
            phi_event = ( phi_low + MCvect[1]*(phi_hi-phi_low) );
        }
        else{
            theta_lab = TMath::ACos(rand->Uniform(cos_theta_hi,cos_theta_low));
            theta_event = TMath::Pi() - theta_lab;
            phi_event = rand->Uniform(-1.*TMath::Pi(),TMath::Pi());

            MCvect[0] = (theta_lab - theta_low )/(theta_hi-theta_low);
            MCvect[1] = ( phi_event + TMath::Pi() ) / (2.*TMath::Pi());
            MCwt = diff_cs(num_dim,MCvect);
            MCwt /= TMath::Sin(theta_lab); //Because of foam generation

            cs_sum+=MCwt; //Because of foam generation
        }

        //Calculate energies and momentum
        Eef = ( Mp*Eei ) /( Mp + Eei*(1.-TMath::Cos(theta_lab)) ); //Final Electron Energy
        Ef_rest.SetE(Eef);
        Ef_rest.SetPx( TMath::Sqrt(Eef*Eef - Me*Me)*TMath::Sin(theta_event)*TMath::Cos(phi_event) );
        Ef_rest.SetPy( TMath::Sqrt(Eef*Eef - Me*Me)*TMath::Sin(theta_event)*TMath::Sin(phi_event) );
        Ef_rest.SetPz( TMath::Sqrt(Eef*Eef - Me*Me)*TMath::Cos(theta_event) );

        Pf_rest.SetE( Eei + Mp - Eef );
        Pf_rest.SetPx( -1.*Ef_rest.Px() );
        Pf_rest.SetPy( -1.*Ef_rest.Py() );
        Pf_rest.SetPz( Ei_rest.Pz() - Ef_rest.Pz() );

        //Boost back into lab frame
        Ef = Ef_rest; Ef.Boost(boost_vec);
        Pf = Pf_rest; Pf.Boost(boost_vec);
        q4vec = Ei - Ef;

        //Calculate Q2
        Q2 = 4.*Eei*Eef*TMath::Sin(theta_lab/2.)*TMath::Sin(theta_lab/2.);

        if(iEvent==0){
            cout<<"Event "<<iEvent<<endl;
            cout<<"Q2 = "<<Q2<<" GeV^2"<<endl;
            
            cout<<"Final Electron 4 Momentum in Proton Rest Frame:"<<endl;
            Ef_rest.Print();
            cout<<"Final Proton 4 Momentum in Proton Rest Frame:"<<endl;
            Pf_rest.Print();

            cout<<"Final Electron 4 Momentum in Lab Frame:"<<endl;
            Ef.Print();
            cout<<"Final Proton 4 Momentum in Lab Frame:"<<endl;
            Pf.Print();

            cout<<endl;
        }

        //Write out event information to the HepMC3 file 
        GenParticlePtr p1 = std::make_shared<GenParticle>(
            FourVector(Ei.Px(),Ei.Py(),Ei.Pz(),Ei.E()), 11, 4);
        GenParticlePtr p2 = std::make_shared<GenParticle>(
            FourVector(Pi.Px(),Pi.Py(),Pi.Pz(),Pi.E()), 2212, 4);
        
        GenParticlePtr p3 = std::make_shared<GenParticle>(
            FourVector(Ef.Px(),Ef.Py(),Ef.Pz(),Ef.E()), 11, 1);
        GenParticlePtr p4 = std::make_shared<GenParticle>(
            FourVector(Pf.Px(),Pf.Py(),Pf.Pz(),Pf.E()), 2212, 1);

        GenVertexPtr v1 = std::make_shared<GenVertex>();
        v1->add_particle_in(p1);
        v1->add_particle_in(p2);

        v1->add_particle_out(p3);
        v1->add_particle_out(p4);
        evt.add_vertex(v1);

        //Save the total cross section (in fb) to the HepMC file.    
        std::shared_ptr<GenCrossSection> cross_section = std::make_shared<GenCrossSection>();
        evt.add_attribute("GenCrossSection",cross_section);

        if(do_foam){
            Double_t mcResult,mcError;
            Foam->GetIntegMC(mcResult,mcError);

            auto mcResult_scaled = mcResult*(theta_hi-theta_low)*(phi_hi-phi_low);
            auto mcError_scaled =  mcError*(theta_hi-theta_low)*(phi_hi-phi_low);
            cross_section->set_cross_section(mcResult_scaled,mcError_scaled);
        }
        else{
            cross_section->set_cross_section(CS_Int_Q2,0.);
        }

        //Save event weight
        std::shared_ptr<Attribute> weight = std::make_shared<DoubleAttribute>(MCwt);
        evt.add_attribute("weight",weight);

        //Print out first event in HepMC format
        if (iEvent==0){
            Print::listing(evt);
            cout<<""<<endl;
        }

        hepmc_output.write_event(evt);
        evt.clear();
    }

    //Total Integral
    if(do_foam==1){
        Double_t mcResult,mcError;
        Foam->GetIntegMC(mcResult,mcError);
        cout <<endl<<"mcResult = "<<mcResult<<" +- "<<mcError<<endl;
        cout <<"Total Cross Section (mcResult Scaled) = "<<mcResult*(theta_hi-theta_low)*(phi_hi-phi_low)<<" +- "<<
            mcError*(theta_hi-theta_low)*(phi_hi-phi_low)<<" Femto Barns"<<endl;
        cout <<"Total Cross Section (Q2 Integral) = "<<CS_Int_Q2<<" Femto Barns"<<endl; 

        cout<<"----------------------"<<endl;
        cout<<"Lower Q2 limit = "<<Q2_low<<" GeV^2"<<endl;
        cout<<"Upper Q2 limit = "<<Q2_hi<<" GeV^2"<<endl;
        cout<<"We have generated "<<nevents<<" events. This corresponds to an integrated luminosity of "<<nevents/CS_Int_Q2<<" fb^-1"<<endl;
        cout<<"For 100 fb^-1 intgrated luminosity, generate "<<CS_Int_Q2*100 <<" events"<<endl;
        cout<<"----------------------"<<endl;

    }
    else{
        cout <<endl<<"Total Cross Section (MCResult) = "<<(cs_sum/((Double_t)nevents) )*(phi_hi-phi_low)*(cos_theta_low-cos_theta_hi)<< " Femto Barns"<<endl;
        cout <<"Total Cross Section (Q2 Integral) = "<<CS_Int_Q2<<" Femto Barns"<<endl;
    }

    //Close output HepMC3 file
    hepmc_output.close();

}