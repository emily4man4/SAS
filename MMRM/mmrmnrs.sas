******************************************************************
*
* Program Name   : mmrmnrs.sas
* Protocol/Study : GWPharm ISE
* Type           : Analysis Reporting Dataset
* Decription     : Program to generate analysis dataset containing
*                  MMRM analysis of change from baseline in 
*                  NRS-spasticity score
*                  Analysis is repeated for studies GWMS0106, GWCL0403,
*                  GWCP0604 and SAVANT and further for the two secondary
*                  baseline definitions for SAVANT
*                  See ISE SAP Section 3.5.2
*
* Author         : Emily Foreman
* Date Created   : 11SEP2020
* Input datasets :
* Macros used    :
* Files used     :
*
*********************************************************************************************
* Change History  
*
* Changed by        : Nitish Ramparsad
* Reason for Change : Added missing label attributes as per QC comments emailed on 10AUG2020
* Date changed      : 10AUG2020  
*
* Changed by        : Emily Foreman
* Reason for Change : Updated to only include planned visits
* Date changed      : 28OCT2020 
*
* Changed by        : Emily Foreman
* Reason for Change : Updated to remove MMRMNRS07/08 analyses
* Date changed      : 02NOV2020
*
* Changed by        : Emily Foreman
* Reason for Change : Updated analysis description for MMRM analysis using subset of patients
                      with mean number of sprays <=12 in ALL weeks of treatment
* Date changed      : 29JAN2021
**********************************************************************************************; 

proc datasets library = work kill nolist;
run;
quit;

*Macro to test convergence of model using unstuctured covariance structure. Creates macro variable, conv, that resolves to Y if 
 the convergence criteria is met and N otherwise, in which case the compound symmetry covariance structure should be used;
%macro convergence();

    %*** Check convergence to determine covariance structure ***;
    %*** ODS select statements temporarily suspend any SAS output because this MIXED procedure is only to determine convergence ***;
    ods select none;
    ods output convergencestatus = convstat;
    ods graphics off;
    proc mixed data = work.&label._adss (where = (1<=avisitn<=&final_vis. and chg ne . and trtpn ne .)) method = reml;
        class trtpn usubjid avisitn;
        model chg = base trtpn avisitn trtpn*avisitn base*avisitn / ddfm = KenwardRoger; 
        repeated avisitn / type = un sub = usubjid;
    run;
    ods output close;
    ods select all;

	%*Use convergencestatus dataset to test convergence of the model, status = 0 and pdg = pdh = 1 indicates model has converged;
	data work.convstat1;
	  set work.convstat;
	  if status = 0 and pdg and pdh then convstat = 'Y';
	  else convstat = 'N';
	run;

	%*Create macro variable to indicate model convergence using unstructured covariance structure;
	data _null_;
	  set work.convstat1;
	  call symputx("conv", convstat);
	run;

	%put &conv;

%mend convergence;

*Macro to run the MMRM analysis for each study and baseline definition, i.e. repeat for secondary BL definitions for SAVANT study;
%macro mmrm(analyid = ,
            analydsc = ,
            study = ,
            base  = ,
            basen = 1,
            nmspr12 = N,
			nmspr12p = N,
            final_vis = );

  %if "&study." = "H15/02" %then %do;
    %if "&nmspr12." = "Y" %then %do; %let label = savant_nspr; %end;
	%else %if "&nmspr12p." = "Y" %then %do; %let label = savant_nsprp; %end;
    %else %do; %let label = savant&basen.; %end;
  %end;
  %else %do;
    %if "&nmspr12." = "Y" %then %do; %let label = &study._nspr; %end;
	%else %if "&nmspr12p." = "Y" %then %do; %let label = &study._nsprp; %end;
    %else %do; %let label = &study.; %end;
  %end;
  
  %*Subset data to the NRS-Spasticity Score Weekly Average for the relevant study, population (efficacy analysis set) and baseline definition;
  data work.&label._adss;
    set adam.adss;
    where studyid = "&study." and paramcd = 'NRSPASAV' and basetype = "&base." and efffl = 'Y'
	      %if "&nmspr12." = "Y" %then %do; and spy12yn = 'Y' %end;
		  ;
	%if "&nmspr12p." = "Y" %then %do;
	  if trtpn = 1 then do;
	    if spy12yn = 'Y';
	  end;
    %end;
  run;

  proc sql noprint;
    select *
	from work.&label._adss;
	%let nobs = &sqlobs.;
  quit;

  %if &nobs. ne 0 %then %do;

    %*Calculate the small n counts for each week, in each treatment group. The MMRM model excludes subjects with missing covariates,
      therefore need to also exclude these from the small n count;
    proc freq data = work.&label._adss (where = (1<=avisitn<=&final_vis. and chg ne . and trtpn ne .));
	  tables trtpn*avisitn / out = work.&label._count;
    run;
    
    data work.&label._count_1 (keep = trtn timeptw n);
      length timeptw $ 8;
      set work.&label._count (rename = (trtpn = trtn
	                                    count = n));
	  timeptw = cat("Week ", strip(put(avisitn, 2.0)));
    run;
    
    %*Call the convergence macro as defined above to test convergence using unstructured covariance structure;
    %convergence();
    
    %*If model converges then covariance structure is unstructured
     else if model does not converge then covariance structure is compound symmetric;
    %if &conv = Y %then %do;
    
      ods output lsmeans = work.&label._lsmeans diffs = &label._diffs;
      proc mixed data = work.&label._adss (where = (1<=avisitn<=&final_vis. and chg ne . and trtpn ne .)) method = reml;
        class trtpn avisitn usubjid;
        model chg = base trtpn avisitn trtpn*avisitn base*avisitn / DDFM = KenwardRoger;
        repeated avisitn / type = un subject = usubjid rcorr;
        lsmeans trtpn*avisitn / cl pdiff;
      run;
    
    %end;
    %else %do;
    
      ods output lsmeans = work.&label._lsmeans diffs = &label._diffs;
      proc mixed data = work.&label._adss (where = (1<=avisitn<=&final_vis. and chg ne . and trtpn ne .)) method = reml;
        class trtpn avisitn usubjid;
        model chg = base trtpn avisitn trtpn*avisitn base*avisitn / DDFM = KenwardRoger;
        repeated avisitn / type = cs subject = usubjid rcorr;
        lsmeans trtpn*avisitn / cl pdiff;
      run;
    
    %end;
    
    data work.&label._lsmeans_1 (keep = analyid analydsc studyid timeptw timeptwn trt trtn lsmean lsmeanse lowerci upperci);
      length timeptw $ 8;
      set work.&label._lsmeans (rename = (trtpn = trtn
                                          estimate = lsmean
	  									  stderr = lsmeanse
	  									  lower = lowerci
	  									  upper = upperci));
      analyid = "&analyid.";
	  analydsc = "&analydsc.";
	  studyid = "&study.";
	  timeptwn = avisitn;
	  timeptw = cat("Week ", strip(put(avisitn, 2.0)));
	  if trtn = 1 then trt = 'Nabiximols';
	  else if trtn = 2 then trt = 'Placebo';
    run;
    
    proc sql;
      create table work.&label._lsmeans__merge as
	  select one.*, two.n
	  from work.&label._lsmeans_1 as one left join work.&label._count_1 as two
	  on one.timeptw = two.timeptw and one.trtn = two.trtn;
    quit;
    
    data work.&label._diffs_1 (keep = timeptw trtn lsmdif lsmdlci lsmduci diffpval);
      length timeptw $ 8;
      set work.&label._diffs (rename = (trtpn = trtn
	                                    _trtpn = _trtn
	                                    estimate = lsmdif
	  								    lower = lsmdlci
	  								    upper = lsmduci
	  								    probt = diffpval));
	  timeptw = cat("Week ", strip(put(avisitn, 2.0)));
      where trtn = 1 and _trtn = 2 and avisitn = _avisitn;
    run;
    
    proc sql;
      create table work.&label._merge as
      select one.*, two.lsmdif, two.lsmdlci, two.lsmduci, two.diffpval
	  from work.&label._lsmeans__merge as one left join work.&label._diffs_1 as two
	  on one.timeptw = two.timeptw and one.trtn = two.trtn;
    quit;

  %end;

  %else %do;

    data work.&label._merge;
      analyid = '';
	  analydsc = ' ';
      studyid = '';
      trt = '';
	  trtn = .;
      timeptw = '';
	  timeptwn = .;
	  n = .;
	  lsmean = .;
	  lsmeanse = .;
	  lowerci = .;
	  upperci = .;
	  lsmdif = .;
	  lsmdlci = .;
	  lsmduci = .;
	  diffpval = .;
    run;

  %end;

%mend mmrm;

%mmrm(analyid = mmrmnrs01, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity),                    study = GWMS0106, base = Primary, final_vis = 6)
%mmrm(analyid = mmrmnrs02, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity),                    study = GWCL0403, base = Primary, final_vis = 14)
%mmrm(analyid = mmrmnrs03, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity),                    study = GWSP0604, base = Primary, final_vis = 12)
%mmrm(analyid = mmrmnrs04, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity),                    study = H15/02,   base = Primary, basen = 1, final_vis = 12)
%mmrm(analyid = mmrmnrs05, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity (7 days Baseline)),  study = H15/02,   base = Secondary NRS-spasticity previous 7 days, basen = 2, final_vis = 12)
%mmrm(analyid = mmrmnrs06, analydsc = %str(MMRM Analysis of CFB in the Mean NRS-Spasticity (14 days Baseline)), study = H15/02,   base = Secondary NRS-spasticity previous 14 days, basen = 3, final_vis = 12)
/*%mmrm(analyid = mmrmnrs07, analydsc = %str(MMRM Analysis of CFB in NRS-Spasticity Score in a Subset of patients with Mean Number of Sprays <=12 in the Last Week of Treatment in Nabiximols Group), study = GWMS0106, base = Primary, nmspr12 = Y, final_vis = 6)*/
/*%mmrm(analyid = mmrmnrs08, analydsc = %str(MMRM Analysis of CFB in NRS-Spasticity Score in a Subset of patients with Mean Number of Sprays <=12 in the Last Week of Treatment in Nabiximols Group), study = GWCL0403, base = Primary, nmspr12 = Y, final_vis = 14)*/
%mmrm(analyid = mmrmnrs09, analydsc = %str(MMRM Analysis of CFB in NRS-Spasticity Score in a Subset of patients with Mean Number of Sprays <=12 in All Weeks of Treatment in Nabiximols Group Versus All Placebo Patients), study = GWMS0106, base = Primary, nmspr12p = Y, final_vis = 6)
%mmrm(analyid = mmrmnrs10, analydsc = %str(MMRM Analysis of CFB in NRS-Spasticity Score in a Subset of patients with Mean Number of Sprays <=12 in All Weeks of Treatment in Nabiximols Group Versus All Placebo Patients), study = GWCL0403, base = Primary, nmspr12p = Y, final_vis = 14)

%*Merge the data for each study and baseline definition into one output dataset;
data work.mmrmnrs;
attrib
    analyid   label = "Analysis ID"
    analydsc  label = "Analysis Description"
	pop       label = "Analysis Set"
	studyid   label = "Study Identifier"
	trt       label = "Planned Treatment"
	trtn      label = "Planned Treatment (N)"
	timeptw   label = "Analysis Timepoint"
	timeptwn  label = "Analysis Timepoint (N)"
	n         label = "n"
	lsmean    label = "LS Mean"
	lsmeanse  label = "SE of LS Mean"
	lowerci   label = "LS Mean Low CI"
	upperci   label = "LS Mean Upper CI"
	lsmdif    label = "LS Mean Difference"
	lsmdlci   label = "Lower 95% CI of LS Mean Difference"
	lsmduci   label = "Upper 95% CI of LS Mean Difference"
	diffpval  label = "LS Mean Difference P-Value";
  length analyid analydsc pop studyid trt $ 200;
  format trtn best8.;
  length timeptw $ 200;
  format timeptwn n lsmean lsmeanse lowerci upperci lsmdif lsmdlci lsmduci diffpval best8.;
  set work.gwms0106_merge
      work.gwcl0403_merge
	  work.gwsp0604_merge
	  work.savant1_merge
	  work.savant2_merge
	  work.savant3_merge
/*	  work.gwcl0403_nspr_merge*/
/*	  work.gwms0106_nspr_merge*/
	  work.gwcl0403_nsprp_merge
	  work.gwms0106_nsprp_merge;
  pop = 'EFF';
run;

proc sort data = work.mmrmnrs
           out = output.mmrmnrs (label='MMRM CFB NRS Spasticity');
  by analyid studyid timeptwn descending trtn;
run;

%checklog()
