functions{
  vector doField(vector phi, vector theta, real[] x, int[] y){
    
    //integer to define array size
    int numMonths = y[1];
    int numFields = y[2];
    
    //derived parameters
    real FracSolidXXXToHUMv263 = 1.0 - phi[6]; //0.54
    real FracSolidHUMToHUMv263 = 1.0 - phi[7]; //0.54
    real Depth = 25.0; //cm
    
    //Calculations Preliminary
    real FracPlantDPM = phi[14] / (1.0 + phi[14]);
    real FracPlantRPM = 1.0 / (1.0 + phi[14]);
    real RatioCO2ToSolids = 1.67*(1.85 + 1.6*exp(-7.86*phi[13]));
    real DcmpFracCO2 = RatioCO2ToSolids / (1.0 + RatioCO2ToSolids);
    real FracXXXToBIOF = phi[6]/(1.0 + RatioCO2ToSolids);
    real FracXXXToBIOS = 0.0/(1.0 + RatioCO2ToSolids);
    real FracXXXToHUM = FracSolidXXXToHUMv263/(1.0 + RatioCO2ToSolids);
    real FracHUMToBIOF = 0.0/(1.0 + RatioCO2ToSolids);
    real FracHUMToBIOS = phi[7]/(1.0 + RatioCO2ToSolids);
    real FracHUMToHUM = FracSolidHUMToHUMv263/(1.0 + RatioCO2ToSolids);
    real DcmpRateBIOF = phi[10]; //per year
    real DcmpRateBIOS = phi[11];//per year
    
    real DPMDcmpFrac;
    real RPMDcmpFrac;
    real BIOFDcmpFrac;
    real BIOSDcmpFrac;
    real HUMDcmpFrac;
    real LastDPMCM;
    real LastRPMCM;
    real LastBIOFCM;
    real LastBIOSCM;
    real LastHUMCM;
    
    real muD;
    real muR;
    real muF;
    real muS;
    real muH;
    
    real thisField_K_R;
    real thisField_K_H;
    real thisField_K_D;
    real thisField_alpha;
    
    // indices for field specific parameters in theta
    int Dindex;
    int Rindex;
    int Findex;
    int Sindex;
    int Hindex;
    
    int Dindex_init;
    int Rindex_init;
    int Findex_init;
    int Sindex_init;
    int Hindex_init;
    int Iindex_init;
    
    int FracManureDPM_index;
    int FracManureRPM_index;
    int FracManureBIOF_index;
    int FracManureBIOS_index;
    int FracManureHUM_index;
    int FracSolidXXXToBIOFv26_index;
    int FracSolidHUMToBIOSv263_index;
    int DcmpRateDPM_index;
    int DcmpRateRPM_index;
    int DcmpRateBIOFv263_index;
    int DcmpRateBIOSv263_index;
    int DcmpRateHUM_index;
    int FracClay_index;
    int RatioDPMToRPM;
    int EvapFactor;
    int sigma2_D;
    int sigma2_R;
    int sigma2_F;
    int sigma2_S;
    int sigma2_H;
    int sigma2_RPM;
    int sigma2_ROC;
    int sigma2_TOC;
    

    
    //indices for the quantities that are in x
    int rateABC_index;
    int Manure_index;
    int PlantResidue_index;
    int mTOC_index;
    int mROC_index;
    int mPOC_index;
    
    //output variable
    real lp;
    
    Dindex_init = 5*numMonths + 1;
    Rindex_init = 5*numMonths + 2;
    Findex_init = 5*numMonths + 3;
    Sindex_init = 5*numMonths + 4;
    Hindex_init = 5*numMonths + 5;
    Iindex_init = 5*numMonths + 6;
      
    //indices for other quantities that are in phi
    
    lp = 0.0;
    
     lp += normal_lpdf(theta[Dindex_init] | 0.0, 0.1);
     lp += normal_lpdf(theta[Rindex_init] | 0.0, 100.0);
     lp += normal_lpdf(theta[Findex_init] | 0.0, 0.01);
     lp += normal_lpdf(theta[Sindex_init] | 0.0, 0.01);
     lp += normal_lpdf(theta[Hindex_init] | 0.0, 100.0);
     lp += normal_lpdf(theta[Iindex_init] | 0.0, 10.0);
    
    
    for(t in 1:numMonths)
    {
      if(t == 1)
      {
        rateABC_index = 1;
        Manure_index = numMonths + 1;
        PlantResidue_index = 2*numMonths + 1;
        mTOC_index = 3*numMonths + 1;
        mROC_index = 4*numMonths + 1;
        mPOC_index = 5*numMonths + 1;
        
        LastDPMCM = theta[Dindex_init];
        LastRPMCM = theta[Rindex_init];
        LastBIOFCM = theta[Findex_init];
        LastBIOSCM = theta[Sindex_init];
        LastHUMCM = theta[Hindex_init];
        
        Dindex = 1;
        Rindex = numMonths + 1;
        Findex = 2*numMonths + 1;
        Sindex = 3*numMonths + 1;
        Hindex = 4*numMonths + 1;
      }
      else
      {
        rateABC_index = rateABC_index + 1;
        Manure_index = Manure_index + 1;
        PlantResidue_index = PlantResidue_index + 1;
        mTOC_index = mTOC_index + 1;
        mROC_index = mROC_index + 1;
        mPOC_index = mPOC_index + 1;
        
        LastDPMCM = theta[Dindex];
        LastRPMCM = theta[Rindex];
        LastBIOFCM = theta[Findex];
        LastBIOSCM = theta[Sindex];
        LastHUMCM = theta[Hindex];
        
        Dindex = Dindex + 1;
        Rindex = Rindex + 1;
        Findex = Findex + 1;
        Sindex = Sindex + 1;
        Hindex = Hindex + 1;
      }
      
      thisField_K_D = phi[8];
      thisField_K_R = phi[9];
      thisField_K_H = phi[12];
      
      
      //alpha terms are actually log alpha terms (treatments)
      if(x[6*numMonths + 1] == 1)
      {
      	//Fallow
		thisField_alpha = exp(phi[37]);
      }
      if(x[6*numMonths + 1] == 2)
      {
      	//Pasture
      	thisField_alpha = exp(phi[24]);
      }
      if(x[6*numMonths + 1] == 3)
      {
      	//NN0
      	thisField_alpha = exp(phi[25]);
      }
      if(x[6*numMonths + 1] == 4)
      {
      	//NN1
      	thisField_alpha = exp(phi[26]);
      }
      if(x[6*numMonths + 1] == 5)
      {
      	//MM0
      	thisField_alpha = exp(phi[27]);
      }
      if(x[6*numMonths + 1] == 6)
      {
      	//MM1
      	thisField_alpha = exp(phi[28]);
      }
      if(x[6*numMonths + 1] == 7)
      {
      	//II0
      	thisField_alpha = exp(phi[29]);
      }
      if(x[6*numMonths + 1] == 8)
      {
      	//II1
      	thisField_alpha = exp(phi[30]);
      }
      if(x[6*numMonths + 1] == 9)
      {
      	//IM0
      	thisField_alpha = exp(phi[31]);
      }
      if(x[6*numMonths + 1] == 10)
      {
      	//IM1
      	thisField_alpha = exp(phi[32]);
      }
      if(x[6*numMonths + 1] == 11)
      {
      	//IN0
      	thisField_alpha = exp(phi[33]);
      }
      if(x[6*numMonths + 1] == 12)
      {
      	//IN1
      	thisField_alpha = exp(phi[34]);
      }
      if(x[6*numMonths + 1] == 13)
      {
      	//MN0
      	thisField_alpha = exp(phi[35]);
      }
      if(x[6*numMonths + 1] == 14)
      {
      	//MN1
      	thisField_alpha = exp(phi[36]);
      }
               
      DPMDcmpFrac = 1.0 - exp(-1*x[rateABC_index] * thisField_K_D * thisField_alpha / 12.0);
      RPMDcmpFrac = 1.0 - exp(-1*x[rateABC_index] * thisField_K_R * thisField_alpha / 12.0);
      BIOFDcmpFrac = 1.0 - exp(-1*x[rateABC_index] * DcmpRateBIOF * thisField_alpha/ 12.0);
      BIOSDcmpFrac = 1.0 - exp(-1*x[rateABC_index] * DcmpRateBIOS * thisField_alpha / 12.0);
      HUMDcmpFrac = 1.0 - exp(-1*x[rateABC_index] * thisField_K_H * thisField_alpha / 12.0);
      
      
      muD = log((1.0 - DPMDcmpFrac) * LastDPMCM + FracPlantDPM * x[PlantResidue_index] + phi[1] * x[Manure_index]) - 0.5*phi[16];
      muR = log((1.0 - RPMDcmpFrac) * LastRPMCM + FracPlantRPM * x[PlantResidue_index] + phi[2] * x[Manure_index])- 0.5*phi[17];
      muF = log((1.0 - BIOFDcmpFrac) * LastBIOFCM + FracXXXToBIOF * (DPMDcmpFrac * LastDPMCM + RPMDcmpFrac * LastRPMCM + BIOFDcmpFrac * LastBIOFCM + BIOSDcmpFrac * LastBIOSCM) + FracHUMToBIOF * HUMDcmpFrac * LastHUMCM + phi[3] * x[Manure_index]) - 0.5*phi[18];
      muS = log((1.0 - BIOSDcmpFrac) * LastBIOSCM + FracXXXToBIOS * (DPMDcmpFrac * LastDPMCM + RPMDcmpFrac * LastRPMCM + BIOFDcmpFrac * LastBIOFCM + BIOSDcmpFrac * LastBIOSCM) + FracHUMToBIOS * HUMDcmpFrac * LastHUMCM + phi[4] * x[Manure_index]) - 0.5*phi[19];
      muH = log((1.0 - HUMDcmpFrac) * LastHUMCM + FracXXXToHUM * (DPMDcmpFrac * LastDPMCM + RPMDcmpFrac * LastRPMCM + BIOFDcmpFrac * LastBIOFCM +BIOSDcmpFrac * LastBIOSCM) + FracHUMToHUM * HUMDcmpFrac * LastHUMCM + phi[5] * x[Manure_index]) - 0.5*phi[20];
      
      
      lp += lognormal_lpdf( theta[Dindex] | muD, sqrt(phi[16]));
      lp += lognormal_lpdf( theta[Rindex] | muR, sqrt(phi[17]));
      lp += lognormal_lpdf( theta[Findex] | muF, sqrt(phi[18]));
      lp += lognormal_lpdf( theta[Sindex] | muS, sqrt(phi[19]));
      lp += lognormal_lpdf( theta[Hindex] | muH, sqrt(phi[20]));
      
      if(x[mTOC_index] != -9999)
      {
        lp += lognormal_lpdf(x[mTOC_index] | log(theta[Dindex] + theta[Rindex] + theta[Findex] + theta[Sindex] + theta[Hindex] + theta[Iindex_init]) - 0.5*phi[23], sqrt(phi[23]));
      }
      
      if(x[mROC_index] != -9999)
      {
        lp += lognormal_lpdf(x[mROC_index] | log(theta[Iindex_init]) - 0.5*phi[22], sqrt(phi[22]));
      }
      
      if(x[mPOC_index] != -9999)
      {
        lp += lognormal_lpdf(x[mPOC_index] | log(theta[Rindex] + theta[Dindex] + theta[Findex]) - 0.5*phi[21], sqrt(phi[21]));
      }
    }
    
    return[lp]';
  } //end doField
} //end functions


data {
    int<lower = 0> numMonths;
    int<lower = 0> numFields;
    matrix<lower = 0>[numFields, numMonths] PlantResidue;
    matrix<lower = 0>[numFields, numMonths] Manure;
    matrix<lower = 0>[numFields, numMonths] RateABC;
    matrix[numFields, numMonths] mTOC;
    matrix[numFields, numMonths] mPOC; //POC/RPM
    matrix[numFields, numMonths] mROC;
    matrix[numFields, 1] treatments;
    //matrix<lower = 0>[numFields, numTimes] mHOC;
}


transformed data {
    //indices 1:numMonths are for rateABC values
    //indices (numMonths + 1):(2*numMonths) are for Manure inputs
    //indices (2*numMonths + 1):(3*numMonths) are for Plant inputs
    //indices (3*numMonths + 1):(4*numMonths) are for mTOC
    //indices (4*numMonths + 1):(5*numMonths) are for mROC
    //indices (5*numMonths + 1):(6*numMonths) are for mPOC
    real shards_x[numFields, 6*numMonths + 3];
    int shards_y[numFields, 2];
    for(f in 1:numFields)
    {
      shards_y[f, 1] = numMonths;
      shards_y[f, 2] = numFields;
      for(i in 1:numMonths)
      {
        shards_x[f, i] = RateABC[f, i];
      }
      for(i in (numMonths + 1):(2*numMonths))
      {
        shards_x[f, i] = Manure[f, (i - numMonths)];
      }
      for(i in (2*numMonths + 1):(3*numMonths))
      {
        shards_x[f, i] = PlantResidue[f, (i - 2*numMonths)];
      }
      for(i in (3*numMonths + 1):(4*numMonths))
      {
        shards_x[f, i] = mTOC[f, (i - 3*numMonths)];
      }
      for(i in (4*numMonths + 1):(5*numMonths))
      {
        shards_x[f, i] = mROC[f, (i - 4*numMonths)];
      }
      for(i in (5*numMonths + 1):(6*numMonths))
      {
        shards_x[f, i] = mPOC[f, (i - 5*numMonths)];
      }
      shards_x[f, (6*numMonths) + 1] = treatments[f, 1];
    }
}



parameters {
  
      //Constants
      real<lower = 0.0, upper = 1.0> FracManureDPM;
      real<lower = 0.0, upper = 1.0> FracManureRPM;
      real<lower = 0.0, upper = 1.0> FracManureBIOF;
      real<lower = 0.0, upper = 1.0> FracManureBIOS;
      real<lower = 0.0, upper = 1.0> FracManureHUM;
      real<lower = 0.0, upper = 1.0> FracSolidXXXToBIOFv263;
      
      real<lower = 0.0, upper = 1.0> FracSolidHUMToBIOSv263;
      
      real<lower = 5.0, upper = 20.0> DcmpRateDPM; //per year
      real<lower = 0.05, upper = 5.0> DcmpRateRPM; //per year
      real<lower = 0.3, upper = 1.0> DcmpRateBIOFv263; //per year
      real<lower = 0.3, upper = 1.0> DcmpRateBIOSv263; //per year
      real<lower = 0.005, upper = 0.05> DcmpRateHUM; //per year
  
      //General
      //Version = 26.3 RothC
      real<lower = 0.0, upper = 1.0> FracClay;
      real<lower = 1.0, upper = 5.0> RatioDPMToRPM;
      real<lower = 0.0, upper = 1.0> EvapFactor;

      matrix<lower = 0, upper = 200>[numFields, numMonths] D;
      matrix<lower = 0, upper = 200>[numFields, numMonths] R;
      matrix<lower = 0, upper = 200>[numFields, numMonths] F;
      matrix<lower = 0, upper = 200>[numFields, numMonths] S;
      matrix<lower = 0, upper = 200>[numFields, numMonths] H;
  
      vector<lower = 0, upper = 200>[numFields] DPMCMInit;
      vector<lower = 0, upper = 200>[numFields] RPMCMInit;
      vector<lower = 0, upper = 200>[numFields] BIOFCMInit;
      vector<lower = 0, upper = 200>[numFields] BIOSCMInit;
      vector<lower = 0, upper = 200>[numFields] HUMCMInit;
      vector<lower = 0, upper = 200>[numFields] IOMCMInit;
      
      real<lower = 0.0, upper = 2.0> sigma2_D;
      real<lower = 0.0, upper = 2.0> sigma2_R;
      real<lower = 0.0, upper = 2.0> sigma2_F;
      real<lower = 0.0, upper = 2.0> sigma2_S;
      real<lower = 0.0, upper = 2.0> sigma2_H;
      
      real<lower = 0.0, upper = 0.5> sigma2_RPM;
      real<lower = 0.0, upper = 0.5> sigma2_ROC;
      real<lower = 0.0, upper = 0.5> sigma2_TOC;
      
      real<lower = -5.0, upper = 5.0> treatment_alpha_Pasture;
      real<lower = -5.0, upper = 5.0> treatment_alpha_Fallow;
      real<lower = -5.0, upper = 5.0> treatment_alpha_NN0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_NN1;
	  real<lower = -5.0, upper = 5.0> treatment_alpha_MM0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_MM1;
      real<lower = -5.0, upper = 5.0> treatment_alpha_II0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_II1;
      real<lower = -5.0, upper = 5.0> treatment_alpha_IM0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_IM1;
      real<lower = -5.0, upper = 5.0> treatment_alpha_IN0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_IN1;
      real<lower = -5.0, upper = 5.0> treatment_alpha_MN0;
      real<lower = -5.0, upper = 5.0> treatment_alpha_MN1;
}


transformed parameters {
  // theta contains the parameters specific to each field
  vector[5*numMonths + 6] shards_theta[numFields];
  // phi contains the global parameters
  vector[37] shards_phi;
  
  // ----- PHI -----
    shards_phi[1] = FracManureDPM / (FracManureDPM + FracManureRPM + FracManureBIOF + FracManureBIOS + FracManureHUM);
    shards_phi[2] = FracManureRPM / (FracManureDPM + FracManureRPM + FracManureBIOF + FracManureBIOS + FracManureHUM);
    shards_phi[3] = FracManureBIOF / (FracManureDPM + FracManureRPM + FracManureBIOF + FracManureBIOS + FracManureHUM);
    shards_phi[4] = FracManureBIOS / (FracManureDPM + FracManureRPM + FracManureBIOF + FracManureBIOS + FracManureHUM);
    shards_phi[5] = FracManureHUM / (FracManureDPM + FracManureRPM + FracManureBIOF + FracManureBIOS + FracManureHUM);
    shards_phi[6] = FracSolidXXXToBIOFv263;
    shards_phi[7] = FracSolidHUMToBIOSv263;
    shards_phi[8] = DcmpRateDPM;
    shards_phi[9] = DcmpRateRPM;
    shards_phi[10] = DcmpRateBIOFv263;
    shards_phi[11] = DcmpRateBIOSv263;
    shards_phi[12] = DcmpRateHUM;
    shards_phi[13] = FracClay;
    shards_phi[14] = RatioDPMToRPM;
    shards_phi[15] = EvapFactor;
    shards_phi[16] = sigma2_D;
    shards_phi[17] = sigma2_R;
    shards_phi[18] = sigma2_F;
    shards_phi[19] = sigma2_S;
    shards_phi[20] = sigma2_H;
    shards_phi[21] = sigma2_RPM;
    shards_phi[22] = sigma2_ROC;
    shards_phi[23] = sigma2_TOC;

    
    shards_phi[24] = treatment_alpha_Pasture;
    shards_phi[25] = treatment_alpha_NN0;
    shards_phi[26] = treatment_alpha_NN1;
	shards_phi[27] = treatment_alpha_MM0;
    shards_phi[28] = treatment_alpha_MM1;
    shards_phi[29] = treatment_alpha_II0;
    shards_phi[30] = treatment_alpha_II1;
    shards_phi[31] = treatment_alpha_IM0;
    shards_phi[32] = treatment_alpha_IM1;
    shards_phi[33] = treatment_alpha_IN0;
    shards_phi[34] = treatment_alpha_IN1;
    shards_phi[35] = treatment_alpha_MN0;
    shards_phi[36] = treatment_alpha_MN1;
    shards_phi[37] = treatment_alpha_Fallow;
    
  
  // ----- THETA -----
  // indices 1:numMonths are the D pool
  // indices (numMonths + 1):(2*numMonths) are the R pool
  // indices (2*numMonths + 1):(3*numMonths) are the F pool
  // indices (3*numMonths + 1):(4*numMonths) are the S pool
  // indices (4*numMonths + 1):(5*numMonths) are for the H pool
  // index (5*numMonths) + 1 is for D_init
  // index (5*numMonths) + 2 is for R_init
  // index (5*numMonths) + 3 is for F_init
  // index (5*numMonths) + 4 is for S_init
  // index (5*numMonths) + 5 is for H_init
  // index (5*numMonths) + 6 is for I_init
  for(f in 1:numFields)
  {
    for(i in 1:numMonths)
    {
      shards_theta[f, i] = D[f, i];
    }
    for(i in (numMonths + 1):(2*numMonths))
    {
      shards_theta[f, i] = R[f, (i - numMonths)];
    }
    for(i in (2*numMonths + 1):(3*numMonths))
    {
      shards_theta[f, i] = F[f, (i - 2*numMonths)];
    }
    for(i in (3*numMonths + 1):(4*numMonths))
    {
      shards_theta[f, i] = S[f, (i - 3*numMonths)];
    }
    for(i in (4*numMonths + 1):(5*numMonths))
    {
      shards_theta[f, i] = H[f, (i - 4*numMonths)];
    }
    shards_theta[f, 5*numMonths + 1] = DPMCMInit[f];
    shards_theta[f, 5*numMonths + 2] = RPMCMInit[f];
    shards_theta[f, 5*numMonths + 3] = BIOFCMInit[f];
    shards_theta[f, 5*numMonths + 4] = BIOSCMInit[f];
    shards_theta[f, 5*numMonths + 5] = HUMCMInit[f];
    shards_theta[f, 5*numMonths + 6] = IOMCMInit[f];
  }
}

model {
      
      //Constants
      FracManureDPM ~ normal(0.49, 0.01);
      FracManureRPM ~ normal(0.49, 0.01);
      FracManureBIOF ~ normal(0.0, 0.01);
      FracManureBIOS ~ normal(0.0, 0.01);
      FracManureHUM ~ normal(0.02, 0.01);
      FracSolidXXXToBIOFv263 ~ normal(0.46, 0.01);
      FracSolidHUMToBIOSv263 ~ normal(0.46, 0.01);
  
      //General
      //Version = 26.3 RothC
      FracClay ~ normal(0.16, 0.02);
      RatioDPMToRPM ~ normal(1.44, 0.5);
      
      EvapFactor ~ normal(0.75, 0.075);
      
      sigma2_D ~ inv_gamma(403.42, 0.3176);
      sigma2_R ~ inv_gamma(403.42, 0.3176);
      sigma2_F ~ inv_gamma(403.42, 0.3176);
      sigma2_S ~ inv_gamma(403.42, 0.3176);
      sigma2_H ~ inv_gamma(403.42, 0.3176);

      sigma2_RPM ~ inv_gamma(10.5, 0.03939044);
      sigma2_ROC ~ inv_gamma(10.5, 0.2899684);
      sigma2_TOC ~ inv_gamma(10.5, 0.05292064);
      
      DcmpRateDPM ~ normal(10.0, 0.5); 
      DcmpRateRPM ~ normal(0.07, 0.0035);
      DcmpRateHUM ~ normal(0.02, 0.001);
      DcmpRateBIOFv263 ~ normal(0.66, 0.033); //per year
      DcmpRateBIOSv263 ~ normal(0.66, 0.033); //per year

	  // these are the log of the multiplicative effects
	  // the fallow treatment is equal to -1 times the sum of the parameters below
      treatment_alpha_Fallow ~ normal(0.0, 1.0);
      treatment_alpha_Pasture ~ normal(0.0, 1.0);
      treatment_alpha_NN0 ~ normal(0.0, 1.0);
      treatment_alpha_NN1 ~ normal(0.0, 1.0);
      treatment_alpha_MM0 ~ normal(0.0, 1.0);
      treatment_alpha_MM1 ~ normal(0.0, 1.0);
      treatment_alpha_II0 ~ normal(0.0, 1.0);
      treatment_alpha_II1 ~ normal(0.0, 1.0);
      treatment_alpha_IM0 ~ normal(0.0, 1.0);
      treatment_alpha_IM1 ~ normal(0.0, 1.0);
      treatment_alpha_IN0 ~ normal(0.0, 1.0);
      treatment_alpha_IN1 ~ normal(0.0, 1.0);
      treatment_alpha_MN0 ~ normal(0.0, 1.0);
      treatment_alpha_MN1 ~ normal(0.0, 1.0);
  
      target += sum(map_rect(doField, shards_phi, shards_theta, shards_x, shards_y));
}

generated quantities {

}
