#include <cmath>
#include <map>
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "Math/GenVector/PxPyPzM4D.h"
#include "Math/GenVector/Boost.h"
#include "TLorentzVector.h"
#include "TH2Poly.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TFile.h"
#include "PhysicsTools/Heppy/interface/Davismt2.h"
#include "TSystem.h"

#include <assert.h>
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"


TString CMSSW_BASE = gSystem->ExpandPathName("${CMSSW_BASE}");

//// UTILITY FUNCTIONS NOT IN TFORMULA ALREADY
//
//
//
//

TF1 * zptFunc = NULL;
TF1 * bdtcdfinvFunc = NULL;

float ptFlavor(float pt, int jetflavor, int targetflavor){

    int ajetflavor = abs(jetflavor);
    if (ajetflavor != targetflavor) return -1.;
    else return pt;

    
}

int convertFlavor(int flavor, float pt =0., float ptcut = 30.){

    int aflav = abs(flavor);
    if (pt > 0. and pt < ptcut) return -2;
    if      (aflav  > 21) return -1;
    else if (aflav == 21) return 0;
    else { return aflav;}

    return -1;
}


int flavBin(int flavor,int nbjet, int min_nbjet){
  if(nbjet<min_nbjet) return -1;
  int aflav = abs(flavor);
  if      (aflav==5) return 2;
  else if (aflav==0) return 0;
  return 1;
}

int ttChannel(bool lep_matched1, bool lep_taufeeddown1, bool lep_matched2, bool lep_taufeeddown2) {
    int code(-1);
    if( (lep_matched1 && !lep_taufeeddown1) && (lep_matched2 && !lep_taufeeddown2) ) code=0;
    if( (lep_matched1 &&  lep_taufeeddown1) && (lep_matched2 && !lep_taufeeddown2) ) code=1;
    if( (lep_matched1 && !lep_taufeeddown1) && (lep_matched2 &&  lep_taufeeddown2) ) code=1;
    if( (lep_matched1 &&  lep_taufeeddown1) && (lep_matched2 &&  lep_taufeeddown2) ) code=2;
    //if( code < 0) std::cout << " this is not leptonic. code -1 for ttbar " << std::endl;
    return code;
}

float sfVariation(float sf1, float sf1unc, int pdg1, float sf2, float sf2unc, int pdg2, int pdg, float var){
    // var should be either 1 or -1
    float newSF = 1.;
    if ( pdg == abs(pdg1) ) newSF = newSF * (sf1 + var*sf1unc) / sf1;
    if ( pdg == abs(pdg2) ) newSF = newSF * (sf2 + var*sf2unc) / sf2;

    // std::cout << "var: " << var << " this is the flavor " << pdg1*pdg2 << " this is the sf1 " << sf1 << " +- " << sf1unc << " this is the sf2 " << sf2 << " +- " << sf2unc << std::endl;
    // std::cout << " this flavor should be non unity: " << pdg << std::endl;
    // std::cout << " new SF before return: " << newSF << std::endl;
    
    return newSF;
}

float weightW(int id1, int id2){
    if      (abs(id1*id2) == 121) return ( 3.58-0.00)/ 945.66;
    else if (abs(id1*id2) == 143) return (10.75-0.82)/7240.57;
    else if (abs(id1*id2) == 169) return ( 0.28-0.00)/2014.89;
    return -999.;
}

float weightMixed(int id1, int id2, int mixrank, int nb = -1){
    float scalefactor = 1.;
    if      (abs(id1*id2) == 121) {
        // blinded if (mixrank == 0) return (3.31-0.12)/23165.;
        // blinded if (mixrank == 1) return (3.31-0.12)/25126.;
        // blinded if (mixrank == 2) return (3.31-0.12)/25836.;
        // blinded if (mixrank == 3) return (3.31-0.12)/26604.;
        // blinded if (mixrank == 4) return (3.31-0.12)/26253.;
        if (mixrank == 0) scalefactor = (19.-0.56)/32478.;
        if (mixrank == 4) scalefactor = (19.-0.56)/35367.;
    }
    else if (abs(id1*id2) == 143) {
        // blinded if (mixrank == 0) return (7.72-0.89)/106.;
        // blinded if (mixrank == 1) return (7.72-0.89)/136.;
        // blinded if (mixrank == 2) return (7.72-0.89)/138.;
        // blinded if (mixrank == 3) return (7.72-0.89)/139.;
        // blinded if (mixrank == 4) return (7.72-0.89)/97.;
        if (mixrank == 0) scalefactor = (26.-3.20)/495.;
        if (mixrank == 4) scalefactor = (26.-3.20)/487.;
    }
    else if (abs(id1*id2) == 169) {
        // blinded if (mixrank == 0) return (3.31-0.78)/72709.;
        // blinded if (mixrank == 1) return (3.31-0.78)/77721.;
        // blinded if (mixrank == 2) return (3.31-0.78)/79937.;
        // blinded if (mixrank == 3) return (3.31-0.78)/81063.;
        // blinded if (mixrank == 4) return (3.31-0.78)/79988.;
        if (mixrank == 0) scalefactor = (11.-2.53)/76531.;
        if (mixrank == 4) scalefactor = (11.-2.53)/82481.;
    }

    if (nb > -0.5){
        if (nb == 0) scalefactor *= 1.0671;
        if (nb == 1) scalefactor *= 0.8558;
        if (nb == 2) scalefactor *= 0.0000;
    }

    return scalefactor;
}

float systWgt(int nb, int syst, int var){

// sorry for the ugly code...

    if        (syst == 1) { // b-tagging b
        if (nb == 0) {
            if (var > 0) return 0.978;
            if (var < 0) return 1.022;
        }
        if (nb == 1) {
            if (var > 0) return 1.011;
            if (var < 0) return 0.990;
        }
        if (nb == 2) {
            if (var > 0) return 1.031;
            if (var < 0) return 0.970;
        }
    } else if (syst == 2) { // b-tagging light
        if (nb == 0) {
            if (var > 0) return 0.999;
            if (var < 0) return 1.003;
        }
        if (nb == 1) {
            if (var > 0) return 1.000;
            if (var < 0) return 1.000;
        }
        if (nb == 2) {
            if (var > 0) return 1.003;
            if (var < 0) return 0.992;
        }
    } else if (syst == 3) { //quenching
        if (nb == 0) {
            if (var > 0) return 1.066;
            if (var < 0) return 1.;
        }
        if (nb == 1) {
            if (var > 0) return 0.977;
            if (var < 0) return 1.;
        }
        if (nb == 2) {
            if (var > 0) return 0.893;
            if (var < 0) return 1.;
        }
    } else if (syst == 4) { // JEC
        if (nb == 0) {
            if (var > 0) return 1.007;
            if (var < 0) return 0.993;
        }
        if (nb == 1) {
            if (var > 0) return 0.998;
            if (var < 0) return 1.001;
        }
        if (nb == 2) {
            if (var > 0) return 0.988;
            if (var < 0) return 1.014;
        }
    } else if (syst == 5) { // JER
        if (nb == 0) {
            if (var > 0) return 1.021;
            if (var < 0) return 0.978;
        }
        if (nb == 1) {
            if (var > 0) return 0.991;
            if (var < 0) return 1.009;
        }
        if (nb == 2) {
            if (var > 0) return 0.965;
            if (var < 0) return 1.038;
        }
    }

    //std::cout << "something went wrong in the systWgtFunction!!! nb: " << nb << " syst: " << syst << " var: " << var << std::endl;
    return 1.;

}

float zNbWeight(int flavor, int nb){

    float weight = 1.;

    if        (flavor == 121){
        if (nb == 0) weight = 0.976;
        if (nb == 1) weight = 1.296;
        if (nb == 2) weight = 1.157;
    } else if (flavor == 143){
        if (nb == 0) weight = 0.957;
        if (nb == 1) weight = 1.463;
        if (nb == 2) weight = 1.178;
    } else if (flavor == 169){
        if (nb == 1) weight = 1.154;
        if (nb == 0) weight = 0.978;
        if (nb == 2) weight = 1.409;
    }

    return weight;
    
}

float zptWeight(float zpt, int flavor, int updn=0){
    // old one with pp MC if (updn == 0) return 1.;
    // old one with pp MC if (!zptFunc) zptFunc = new TF1("zptFunc","1.14*TMath::Erf((x+4.04)/9.83)",0.,1000.);
    // old one with pp MC return zptFunc->Eval(zpt);
    //
    if (flavor == 143) return 1.;

    // old. possibly from pp if (!zptFunc) zptFunc = new TF1("zptFunc","0.703302 + 0.0704182*x - 0.00294974*TMath::Power(x,2) + 0.0000307175*TMath::Power(x,3)",0.,1000.);
    // old. possibly from pp float weight = zpt >= 45. ? zptFunc->Eval(45.) : zptFunc->Eval(zpt);


    // pol3 if (!zptFunc) zptFunc = new TF1("zptFunc","0.655793 + 0.0647376*x - 0.00234078*TMath::Power(x,2) + (2.14129e-5)*TMath::Power(x,3)",0.,1000.);
    if (!zptFunc) zptFunc = new TF1("zptFunc","0.572243 + 0.0888349*x - 0.00406509*TMath::Power(x,2) + (6.32276e-5)*TMath::Power(x,3) - (3.18089e-07)*TMath::Power(x,4)",0.,1000.);

    float weight = zpt >= 90. ? zptFunc->Eval(90.) : zptFunc->Eval(zpt);

    if        (updn ==  0) {
        return weight;
    } else if (updn ==  1) {
        return weight*TMath::Sqrt(weight);
    } else if (updn == -1) {
        return TMath::Sqrt(weight);
    }
    return 1.;
    
}

float bdtcdfinv(float bdt) {

  if(!bdtcdfinvFunc) bdtcdfinvFunc=new TF1("bdtcdfinv",
                                           "0.364366+0.356660*x+0.091563*pow(x,2)+0.089546*pow(x,3)-0.074305*pow(x,4)+0.066103*pow(x,5)+0.128357*pow(x,6)",-1,1);
  return max(min(1.,bdtcdfinvFunc->Eval(bdt)),0.);
}



float iso02(int pdgId, float pt, float iso02, float lep_rho){
    float ue(-1000.);

    float rho0(0.), a(0.), b(0.), c(0.), d(0.), e(0.), f(0.);

    if (pdgId == 11){
        rho0 = 108.5698;
        a =  0.000707;
        b =  0.063990;
        c =  4.113015;
        d =  0.000000;
        e = 17.075343;
    } else if (pdgId == 13) {
        rho0 = 103.57603;
        a =  0.001974;
        b = -0.022796;
        c =  5.830865;
        d =  0.000000;
        e = 40.140637;
    }

    // isoFormula  = 'x<[0] ? [1]*pow(x,2)+[2]*x+[3] : [4]*pow(log(x),2)+[5]*log(x)+[3]+[1]*pow([0],2)-[4]*pow(log([0]),2)+[2]*[0]-[5]*log([0])'
    float value_true  = a*TMath::Power(lep_rho,2) + b*lep_rho + c;
    float value_false = a*TMath::Power(rho0   ,2) + b*rho0    + c + e*TMath::Log(lep_rho) - e*TMath::Log(rho0);

    ue = lep_rho < rho0 ? value_true : value_false ;

    return (iso02 - ue)/pt;

}

float iso03(int pdgId, float pt, float isofull, float isofull30, float lep_rho){

    float ue(-1000.);
    float iso(-1000.);

    if      (pdgId == 11){
        ue = 0.000817*TMath::Power(lep_rho+14.696,2)+0.201661*(lep_rho+14.696);
        iso  = isofull*pt;
        iso += 0.0011*TMath::Power(lep_rho+142.4,2)-0.14*(lep_rho+142.4);
        iso -= ue;
    }

    else if (pdgId == 13) {
        ue = 0.00102*TMath::Power(lep_rho+12.6255,2)+0.18535*(lep_rho+12.6255);
        iso = isofull30 - ue;
    }

    return iso/pt;
}


TFile * lhinfile = NULL;
TH1F * hist_lh_dphi  = NULL;
TH1F * hist_lh_llpt  = NULL;
TH1F * hist_lh_lleta = NULL;

float simpleLH3d(float dphi, float llpt, float lleta){

    if (!lhinfile || !hist_lh_dphi || !hist_lh_llpt || !hist_lh_lleta){
        lhinfile = new TFile("ALL.root","read");
        hist_lh_dphi  = (TH1F*)(lhinfile->Get("dphi_ttbar"));
        hist_lh_llpt  = (TH1F*)(lhinfile->Get("llpt_ttbar"));
        hist_lh_lleta = (TH1F*)(lhinfile->Get("lleta_ttbar"));
        hist_lh_dphi ->Scale(1./hist_lh_dphi ->Integral());
        hist_lh_llpt ->Scale(1./hist_lh_llpt ->Integral());
        hist_lh_lleta->Scale(1./hist_lh_lleta->Integral());
    }

    int bin_dphi = std::max(1, std::min(hist_lh_dphi ->GetNbinsX(), hist_lh_dphi ->FindBin(dphi)));
    int bin_llpt = std::max(1, std::min(hist_lh_llpt ->GetNbinsX(), hist_lh_llpt ->FindBin(dphi)));
    int bin_lleta= std::max(1, std::min(hist_lh_lleta->GetNbinsX(), hist_lh_lleta->FindBin(dphi)));
    
    float prob_dphi  = max(0.001,hist_lh_dphi ->GetBinContent(bin_dphi ));
    float prob_llpt  = max(0.001,hist_lh_llpt ->GetBinContent(bin_llpt ));
    float prob_lleta = max(0.001,hist_lh_lleta->GetBinContent(bin_lleta));

    //std::cout << "returning " << prob_dphi*prob_llpt*prob_lleta << std::endl;

    return prob_dphi*prob_llpt*prob_lleta*25.*25.*25.;

}


TF1 * pdf_dphi  = NULL;
TF1 * pdf_llpt  = NULL;
TF1 * pdf_lleta = NULL;

float pdfLH3d(float dphi, float llpt, float lleta){

    if (!pdf_dphi || !pdf_llpt || !pdf_lleta){
        pdf_dphi  = new TF1("pdf_dphi" ,"[0]*x+[1]");
        pdf_dphi->SetParameter(0,0.0125184);
        pdf_dphi->SetParameter(1,0.0202382);

        pdf_lleta  = new TF1("pdf_lleta" ,"gaus(0)");
        pdf_lleta->SetParameter(0, 1./TMath::Sqrt(2.*TMath::Pi()*TMath::Power(1.08917e+00,2)));
        pdf_lleta->SetParameter(1, 4.77325e-01);
        pdf_lleta->SetParameter(2, 1.08917e+00);

        pdf_llpt = new TF1("pdf_llpt","gaus(0)");
        pdf_llpt->SetParameter(0, 1./TMath::Sqrt(2.*TMath::Pi()*TMath::Power(3.39914e+01,2)));
        pdf_llpt->SetParameter(1, 5.97579e+01);
        pdf_llpt->SetParameter(2, 3.39914e+01);
    }

    float prob_dphi  = pdf_dphi ->Eval(dphi );
    float prob_llpt  = pdf_llpt ->Eval(llpt );
    float prob_lleta = pdf_lleta->Eval(lleta);

    //std::cout << "returning " << prob_dphi << " " << prob_llpt << " " << prob_lleta << std::endl;

    return prob_dphi*prob_llpt*prob_lleta*1000.*3.;

}

float myratio(float num, float denom) {
  if(denom==0) return 0;
  return num/denom;
}

float deltaPhi(float phi1, float phi2) {
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

float if3(bool cond, float iftrue, float iffalse) {
    return cond ? iftrue : iffalse;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
    float deta = std::abs(eta1-eta2);
    float dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

float pt_2(float pt1, float phi1, float pt2, float phi2) {
    phi2 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2), pt2*std::sin(phi2));
}

float mt_2(float pt1, float phi1, float pt2, float phi2) {
    return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

float rapidity_2(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,0.);
    PtEtaPhiMVector p42(pt2,eta2,phi2,0.);
    return (p41+p42).Rapidity();
}


float mass_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float mt2davis(float pt1, float eta1, float phi1, float pt2, float eta2, float phi2, float met, float metphi){
    // NOTE THAT THIS FUNCTION ASSUMES MASSLESS OBJECTS. NOT ADVISED TO USE WITH HEMISPHERES ETC.
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p1(pt1,eta1,phi1,0.);
    PtEtaPhiMVector p2(pt2,eta2,phi2,0.);
    PtEtaPhiMVector mv(met,0.,metphi,0.);
    double a[] = {p1.M(), p1.Px(), p1.Py()};
    double b[] = {p2.M(), p2.Px(), p2.Py()};
    double c[] = {mv.M(), mv.Px(), mv.Py()};

    heppy::Davismt2 mt2obj;
    mt2obj.set_momenta( a, b, c );
    mt2obj.set_mn( 0. );

    float result = (float) mt2obj.get_mt2();
    return result;
}

float phi_2(float pt1, float phi1, float pt2, float phi2) {
    float px1 = pt1 * std::cos(phi1);
    float py1 = pt1 * std::sin(phi1);
    float px2 = pt2 * std::cos(phi2);
    float py2 = pt2 * std::sin(phi2);
    return std::atan2(py1+py2,px1+px2);
}

float phi_3(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3) {
    float px1 = pt1 * std::cos(phi1);
    float py1 = pt1 * std::sin(phi1);
    float px2 = pt2 * std::cos(phi2);
    float py2 = pt2 * std::sin(phi2);
    float px3 = pt3 * std::cos(phi3);
    float py3 = pt3 * std::sin(phi3);
    return std::atan2(py1+py2+py3,px1+px2+px3);
}

float eta_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Eta();
}

float pt_3(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3) {
    phi2 -= phi1;
    phi3 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3), pt2*std::sin(phi2) + pt3*std::sin(phi3));
}


float mass_3(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    return (p41+p42+p43).M();
}


float pt_4(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3, float pt4, float phi4) {
    phi2 -= phi1;
    phi3 -= phi1;
    phi4 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3) + pt4 * std::cos(phi4), pt2*std::sin(phi2) + pt3*std::sin(phi3) + pt4*std::sin(phi4));
}
 
float mass_4(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float pt4, float eta4, float phi4, float m4) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    PtEtaPhiMVector p44(pt4,eta4,phi4,m4);
    return (p41+p42+p43+p44).M();
}

float mt_llv(float ptl1, float phil1, float ptl2, float phil2, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}

float mt_lllv(float ptl1, float phil1, float ptl2, float phil2, float ptl3, float phil3, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptl3*std::cos(phil3) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptl3*std::sin(phil3) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptl3+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}


float mtw_wz3l(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float mZ1, float met, float metphi) 
{
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)) < 0.01) return mt_2(pt3,phi3,met,metphi);
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt2,phi2,met,metphi);
    if (abs(mZ1 - mass_2(pt2,eta2,phi2,m2,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt1,phi1,met,metphi);
    return 0;
}

float u1_2(float met_pt, float met_phi, float ref_pt, float ref_phi) 
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_px + ref_px;
    return (ux*ref_px + uy*ref_py)/ref_pt;
}
float u2_2(float met_pt, float met_phi, float ref_pt, float ref_phi)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_px + ref_px;
    return (ux*ref_py - uy*ref_px)/ref_pt;
}

// reconstructs a top from lepton, met, b-jet, applying the W mass constraint and taking the smallest neutrino pZ
float mtop_lvb(float ptl, float etal, float phil, float ml, float met, float metphi, float ptb, float etab, float phib, float mb) 
{
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > PxPyPzMVector;
    PtEtaPhiMVector p4l(ptl,etal,phil,ml);
    PtEtaPhiMVector p4b(ptb,etab,phib,mb);
    double MW=80.4;
    double a = (1 - std::pow(p4l.Z()/p4l.E(), 2));
    double ppe    = met * ptl * std::cos(phil - metphi)/p4l.E();
    double brk    = MW*MW / (2*p4l.E()) + ppe;
    double b      = (p4l.Z()/p4l.E()) * brk;
    double c      = met*met - brk*brk;
    double delta   = b*b - a*c;
    double sqdelta = delta > 0 ? std::sqrt(delta) : 0;
    double pz1 = (b + sqdelta)/a, pz2 = (b - sqdelta)/a;
    double pznu = (abs(pz1) <= abs(pz2) ? pz1 : pz2);
    PxPyPzMVector p4v(met*std::cos(metphi),met*std::sin(metphi),pznu,0);
    return (p4l+p4b+p4v).M();
}

float DPhi_CMLep_Zboost(float l_pt, float l_eta, float l_phi, float l_M, float l_other_pt, float l_other_eta, float l_other_phi, float l_other_M){
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
  PtEtaPhiMVector l1(l_pt,l_eta,l_phi,l_M);
  PtEtaPhiMVector l2(l_other_pt,l_other_eta,l_other_phi,l_other_M);
  PtEtaPhiMVector Z = l1+l2;
  ROOT::Math::Boost boost(Z.BoostToCM());
  l1 = boost*l1;
  return deltaPhi(l1.Phi(),Z.Phi());
}


//PU weights

// for json up to 276811 (12.9/fb), pu true reweighting
float _puw2016_nTrueInt_13fb[60] = {0.0004627598152210959, 0.014334910915287028, 0.01754727657726197, 0.03181477917631854, 0.046128282569231016, 0.03929080994013006, 0.057066019809589925, 0.19570744862221007, 0.3720256062526554, 0.6440076202772811, 0.9218024454406528, 1.246743510634073, 1.5292543296414058, 1.6670061646418215, 1.7390553377117133, 1.6114721876895595, 1.4177294439817985, 1.420132866045718, 1.3157656415540477, 1.3365188060918483, 1.1191478126677334, 0.9731079434848392, 0.9219564145009487, 0.8811793391804676, 0.7627315352977334, 0.7265186492688713, 0.558602385324645, 0.4805954159733825, 0.34125298049234554, 0.2584848657646724, 0.1819638766151892, 0.12529545619337035, 0.11065705912071645, 0.08587356267495487, 0.09146322371620583, 0.11885517671051576, 0.1952483711863489, 0.23589115679998116, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw2016_nTrueInt_13fb(int nTrueInt) { if (nTrueInt<60) return _puw2016_nTrueInt_13fb[nTrueInt]; else return 0; }

float _puw2016_nTrueInt_36fb[100] = {0.3505407355600995, 0.8996968628890968, 1.100322319466069, 0.9562526765089195, 1.0366251229154624, 1.0713954619016586, 0.7593488199769544, 0.47490309461978414, 0.7059895997695581, 0.8447022252423783, 0.9169159386164522, 1.0248924033173097, 1.0848877947714115, 1.1350984224561655, 1.1589888429954602, 1.169048420382294, 1.1650383018054549, 1.1507200023444994, 1.1152571438041776, 1.0739529436969637, 1.0458014000030829, 1.032500407707141, 1.0391236062781293, 1.041283620738903, 1.0412963370894526, 1.0558823002770783, 1.073481674823461, 1.0887053272606795, 1.1041701696801014, 1.123218903738397, 1.1157169321377927, 1.1052520327174429, 1.0697489590429388, 1.0144652740600584, 0.9402657069968621, 0.857142825520793, 0.7527112615290031, 0.6420618248685722, 0.5324755829715156, 0.4306470627563325, 0.33289171600176093, 0.24686361729094983, 0.17781595237914027, 0.12404411884835284, 0.08487088505600057, 0.056447805688061216, 0.03540829360547507, 0.022412461576677457, 0.013970541270658443, 0.008587896629717911, 0.004986410514292661, 0.00305102303701641, 0.001832072556146534, 0.0011570757619737708, 0.0008992999249003301, 0.0008241241729452477, 0.0008825716073180279, 0.001187003960081393, 0.0016454104270429153, 0.0022514113879764414, 0.003683196037880878, 0.005456695951503178, 0.006165248770884191, 0.007552675218762607, 0.008525338219226993, 0.008654690499815343, 0.006289068906974821, 0.00652551838513972, 0.005139581024893171, 0.005115751962934923, 0.004182527768384693, 0.004317593022028565, 0.0035749335962533355, 0.003773660372937113, 0.002618732319396435, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw2016_nTrueInt_36fb(int nTrueInt) { if (nTrueInt<100) return _puw2016_nTrueInt_36fb[nTrueInt]; else return 0; }

float mass_3_cheap(float pt1, float eta1, float pt2, float eta2, float phi2, float pt3, float eta3, float phi3) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,0,   0.0);
    PtEtaPhiMVector p42(pt2,eta2,phi2,0.0);
    PtEtaPhiMVector p43(pt3,eta3,phi3,0.0);
    return (p41+p42+p43).M();
}

float lnN1D_p1(float kappa, float x, float xmin, float xmax) {
    return std::pow(kappa,(x-xmin)/(xmax-xmin));
}


void functions() {}





