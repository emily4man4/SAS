****************************************************************************************************
*
* Macro Name        : qc_m_eff207.sas
* Protocol / Study  : D9480C00001 (Harmonize Asia)
* Type              : Macro
* Description       : To validate the Subgroup analysis of the proportion of normokalemic patients at 
*                     Day 29 for Table 14.2.12.2 (eff207), Table 14.2.12.3 (eff207a) Table 14.2.12.4
*                     (eff207b), Table 14.2.12.4 (eff207c)
*
* Author            : Collette Letham
* Date Completed    : 11Feb2022
* 
* Input Dataset     : ADSK
* Macros Used       : %qc_bign, %validation %checklog
*
****************************************************************************************************
* Change History
*
* Amended By      : Emily Foreman
* Date Amended    : 31May2022
* Amended         : Added condition that we require a minimum of 5 responses per treatment-subgroup
*                   for subgroup analysis to be performed
*
* Amended By      : Emily Foreman
* Date Amended    : 11OCT2022
* Amended         : Updated subgroup condition for min number of subjects required to run analysis
*
* Amended By      : Emily Foreman
* Date Amended    : 03NOV2022
* Amended         : Updated to check the number of dps to use when presenting the 95% CI
*
* Amended By      : Emily Foreman
* Date Amended    : 25NOV2022
* Amended         : Updated subgroup condition to use number of events instead of number of subjects
*
****************************************************************************************************;

**Set program name for compare;
%macro qc_m_eff207(subvar    = ,
                   progname  = ,
				   raasi_ref = ,
				   hf_ref    = ,
				   diab_ref  = ,
				   ckd_ref   = );

  %**Calculate Big Ns;
  %qc_bign (pop =FAS02FL);

  **Read ADSK data**;
  data work.norm1;
    set adam.adsk (where=(fas02fl = "Y" and paramcd='SKRESP' and anl01fl = 'Y' and avisitn = 14)) ;

    %*Check for model variables that are missing;
    if cmiss(trt02pn, blgfr1, blsk1, blsk2, agegr1n, braasifl, bckdfl, bhffl, bdiabfl) ne 0 then mis_flag = "Y";
  run;

  %*Check for duplicate records;
  proc sort data = work.norm1 
             out = work.norm1s nodupkey dupout = work.dups;
    by trt02pn trt02p &subvar. usubjid;
  run;

  data _null_;
    set work.dups;
	if _n_ = 1 then put "WAR" "NING: Duplicate records for Day 29 in ADSK, need to investigate why";
  run;

  %*Calculate Big N for each subgroup;
  proc freq data = adam.adsl (where = (fas02fl = 'Y' and trt02pn ne .)) noprint;
    table trt02pn*&subvar. / out = work.bign (drop=percent);
  run;

  %**Small n and response count;
  proc freq data = work.norm1s (where = (mis_flag ne "Y")) noprint;
    table trt02pn*&subvar.       / out = work.n1 (drop=percent);
    table trt02pn*&subvar.*avalc / out = work.n_norm (drop=percent);
  run;

  %*Create Dummy data set to check all the counts have at least 1 observation;
  data work.dummy;
    do trt02pn = 1 to 3;
	  length &subvar. $2;
	  do &subvar. = "N", "Y";
	    length avalc $200;
	    do avalc = "Normokalemic", "Not Normokalemic";
		  output;
		end;
	  end;
	end;
  run;

  %*Merge on counts;
  data work.n_norm2;
    merge work.dummy (in=a) 
          work.n_norm (in=b);
	by trt02pn &subvar. avalc;
	if a and not b then count = 0;
  run;

  %*Merge n counts together for table;
  data work.n_counts;
    merge work.bign (rename=(count=bign)) 
          work.n1 (in=a rename = (count=n1)) 
          work.n_norm2 (in=b rename=(count=n2) where =(avalc = "Normokalemic"));
    by trt02pn &subvar.;
  run;

  data work.n_counts2 (keep = trt02pn &subvar. n2 col: ord1);
    set work.n_counts;
    length col1-col4 $200 id $5 ;

	%*Col1 - Subgroup;
	%if       %upcase(&subvar.) = BCKDFL   %then %do; id = "CKD";   %end;
	%else %if %upcase(&subvar.) = BDIABFL  %then %do; id = "DM";    %end;
	%else %if %upcase(&subvar.) = BHFFL    %then %do; id = "HF";    %end;
	%else %if %upcase(&subvar.) = BRAASIFL %then %do; id = "RAASi"; %end;
	
    if &subvar. = "Y" then do; ord1 = 1; col1 = strip(id); end;
	else if &subvar. = "N" then do; ord1 = 2; col1 = "Not "||strip(id); end;

    %*Col2;
    if trt02pn = 1 then col2 = "ZS 5g (QD)#nN="||strip(put(bign,3.0));
    else if trt02pn = 2 then col2 = "ZS 10g (QD)#nN="||strip(put(bign,3.0));
    else if trt02pn = 3 then col2 = "Placebo#nN="||strip(put(bign,3.0));

    %*Col3;
    if not missing(n1) then col3 = put(n1,2.0);
    else col3 = put(0, 2.0);

    %*Col4;
    if n2 ne 0 then do;
      perc = round(((n2/n1)*100),0.1);
      if perc = 100 then col4 = put(n2, 2.0)||" ( 100)";
	    else col4 = put(n2, 2.0)||" ("||put(perc,4.1)||")";
    end;
    else col4 = put(0, 2.0);
  run;

  %*Only run analysis if there are a minimum of 10 subjects per treatment/subgroup combination;
  data work.check;
    set norm1s (where = (mis_flag ne "Y" and avalc = "Normokalemic"));
	if trt02pn = 1 then do;
	  treat = 1;
	  output;
	end;
	else if trt02pn = 2 then do;
	  treat = 2;
	  output;
	end;
	else if trt02pn = 3 then do;
	  treat = 1;
	  output;
	  treat = 2;
	  output;
	end;
  run;

  proc sql noprint;
    select treat, &subvar., count(*) 
    into: treat1 -,
        : subvar1 -,
        : check1-
    from work.check
    group by treat, &subvar.;
  quit;

  %put &check1. &check2. &check3. &check4.;

  data work.norm2;
    set work.norm1s;
	%*If both treatment combinations for the subgroup level have less than 10 events then remove all records from this subgroup level from the analysis;
	if &check1. <10 and &check3. <10 then do;
	  if &subvar. = "&subvar1." then delete;
    end;
	%*If just one treatment combination for the subgroup level has less than 10 events then remove all records from the corresponding treatment group;
	else do;
	  if &check1. <10 then do;
	    if trt02pn = 1 and &subvar. = "&subvar1." then delete;
	  end;
	  if &check3. <10 then do;
	    if trt02pn = 2 and &subvar. = "&subvar3." then delete;
	  end;
	end;

	%*If both treatment combinations for the subgroup level have less than 10 events then remove all records from this subgroup level from the analysis;
	if &check2. <10 and &check4. <10 then do;
	  if &subvar. = "&subvar2." then delete;
    end;
	%*If just one treatment combination for the subgroup level has less than 10 events then remove all records from the corresponding treatment group;
	else do;
	  if &check2. <10 then do;
	    if trt02pn = 1 and &subvar. = "&subvar2." then delete;
	  end;
	  if &check4. <10 then do;
	    if trt02pn = 2 and &subvar. = "&subvar4." then delete;
	  end;
	end;
  run;

  %*If we have had to remove all records from one level of the subgroup, then ensure the remaining subgroup level is set as the reference otherwise the analysis will not run;
  %if &check1. <10 and &check3. <10 %then %do; %let ref = Y; %end;
  %else %if &check2. <10 and &check4. <10 %then %do; %let ref = N; %end;
  %else %do; %let ref = ; %end;

  %*If there are no treatment combinations in either subgroup level with at least 10 events then do not run the analysis;
  %if &check1.<10 and &check2.<10 and &check3.<10 and &check4.<10  %then %do;

    data work.qc_&progname.;
      length col1 $200;
      col1 = 'NO DATA AVAILABLE FOR THIS REPORT';
    run;

  %end;

  %*If all records from one level of the subgroup has been removed then the following will run;
  %else %if "&ref." = "N" or "&ref." = "Y" %then %do;

	%if       %upcase(&subvar.) = BCKDFL   %then %do; %let ckd_ref = &ref.;   %end;
	%else %if %upcase(&subvar.) = BDIABFL  %then %do; %let diab_ref = &ref.;  %end;
	%else %if %upcase(&subvar.) = BHFFL    %then %do; %let hf_ref = &ref.;    %end;
	%else %if %upcase(&subvar.) = BRAASIFL %then %do; %let raasi_ref = &ref.; %end;

    %*Calculate Odds ratio; 
    proc logistic data = work.norm2 (where = (mis_flag ne "Y"));
      class trt02pn (ref='3') agegr1n braasifl (ref = "&raasi_ref.") bckdfl (ref = "&ckd_ref.") bhffl(ref = "&hf_ref.") bdiabfl(ref = "&diab_ref.") / param=ref;
      model avalc(event = "Normokalemic") = trt02pn blgfr1 blsk1 blsk2 agegr1n braasifl bckdfl bhffl bdiabfl trt02pn*&subvar. / orpvalue;
      oddsratio trt02pn;
      ods output OddsRatiosWald=or1 ParameterEstimates=est1;
    run;

	%*EF 03NOV2022: Added code to check the number of dps required when presenting the 95% CI;
	proc sql noprint;
	  select max(lowercl), max(uppercl)
	  into: max_lower,
          : max_upper
	  from work.or1;
	quit;

	data _null_;
      if &max_lower. ge 1000 then call symputx('lower_dp', 7.2);
      else if &max_lower. ge 100 then call symputx('lower_dp', 6.2);
      else if &max_lower. ge 10 then call symputx('lower_dp', 5.2);
	  else call symputx('lower_dp', 4.2);

      if &max_upper. ge 1000 then call symputx('upper_dp', 7.2);
      else if &max_upper. ge 100 then call symputx('upper_dp', 6.2);
      else if &max_upper. ge 10 then call symputx('upper_dp', 5.2);
	  else call symputx('upper_dp', 4.2);
	run;

    data work.or2;
      length col5-col7 $200;
      set work.or1;

	  %*Set Subvar;
	  if substr(reverse(effect),1,1) = "N" then &subvar. = "N";
	  else if substr(reverse(effect),1,1) = "Y" then &subvar. = "Y";

	  %*Set Treatment and delete rows for ZS 5g v ZS 10g;
	  if substr(effect,1,14) = "TRT02PN 1 vs 2" then delete;
	  else if substr(effect,1,14) = "TRT02PN 1 vs 3" then trt02pn = 1;
	  else if substr(effect,1,14) = "TRT02PN 2 vs 3 " then trt02pn = 2;
  
      if not missing(oddsratioest) and oddsratioest <0.001 then col5 = "<0.001";
      else if not missing(oddsratioest) then col5 = put(oddsratioest, 5.2);

      if nmiss(lowercl, uppercl)=0 then col6 = cat("(", put(lowercl, &lower_dp.), ", ", put(uppercl, &upper_dp.), ")");
      else if missing(lowercl) and not missing(uppercl) then col6 = cat("(   NC, ", put(uppercl, &upper_dp.), ")");
      else if missing(lowercl) and not missing(uppercl) then col6 = cat("(", put(lowercl, &lower_dp.), ", NC)");
	  else if nmiss(lowercl, uppercl)=2 then col6 ="NC";

      if pvalue gt 0.999 then col7 = ">0.999";
      else if pvalue lt 0.001 then col7 = "<0.001";
      else col7 = put(pvalue,6.3);
    run;

    proc sort data = work.or2 out = work.or2s;
      by trt02pn &subvar.;
    run;

    %*Treatment-by-subgroup interaction: p-value only;
    data work.interaction1;
	  length col8 $200;
      set work.est1 (where =(variable=upcase("TRT02PN*&subvar") and classval0 in ("1" "2") and not missing(classval1)));

	  %*Set Treatment;
      if not missing(classval0) then trt02pn = input(classval0,2.);

	  %*Set Subvar;
	  if not missing(classval1) then &subvar = strip(classval1);

	  %*Format p-value;
	  if not missing(probchisq) and probchisq gt 0.999 then col8 = ">0.999";
      else if not missing(probchisq) and probchisq lt 0.001 then col8 = "<0.001";
      else if not missing(probchisq) then col8 = put(probchisq,6.3);
	  else if missing(probchisq) then col8 = "NC";
	
      keep trt02pn &subvar. col8;
    run;

    proc sort data = work.interaction1 out = work.interaction1s;
      by trt02pn &subvar.;
    run;

    %*Merge counts and odds ratio together;
	data work.stats1;
      merge work.n_counts2 
            work.or2s (keep = &subvar. trt02pn col:) 
            work.interaction1s;
      by trt02pn &subvar.;
     
	  if trt02pn ne 3 then do;
	    if missing(col5) then col5 = "  NC";
	    if missing(col6) then col6 = "NC";
	    if missing(col7) then col7 = "NC";
	    if scan(col1, 1, " ") ne "Not" and missing(col8) then col8 = "NC";
	  end;
    run;

    proc sort data = work.stats1 
               out = work.qc_&progname. (keep = col:);
      by ord1 trt02pn;
    run;

  %end;

  %*If there is at least one treatment combination in both subgroup levels with at least 10 events then the following will run;
  %else %do;

    %*Calculate Odds ratio; 
    proc logistic data = work.norm2 (where = (mis_flag ne "Y"));
      class trt02pn (ref='3') agegr1n braasifl (ref = 'N') bckdfl (ref = 'N') bhffl(ref = 'N') bdiabfl(ref = 'N') / param=ref;
      model avalc(event = "Normokalemic") = trt02pn blgfr1 blsk1 blsk2 agegr1n braasifl bckdfl bhffl bdiabfl trt02pn*&subvar. / orpvalue;
      oddsratio trt02pn;
      ods output OddsRatiosWald=or1 ParameterEstimates=est1;
    run;

	%*EF 03NOV2022: Added code to check the number of dps required when presenting the 95% CI;
	proc sql noprint;
	  select max(lowercl), max(uppercl)
	  into: max_lower,
          : max_upper
	  from work.or1;
	quit;

	data _null_;
      if &max_lower. ge 1000 then call symputx('lower_dp', 7.2);
      else if &max_lower. ge 100 then call symputx('lower_dp', 6.2);
      else if &max_lower. ge 10 then call symputx('lower_dp', 5.2);
	  else call symputx('lower_dp', 4.2);

      if &max_upper. ge 1000 then call symputx('upper_dp', 7.2);
      else if &max_upper. ge 100 then call symputx('upper_dp', 6.2);
      else if &max_upper. ge 10 then call symputx('upper_dp', 5.2);
	  else call symputx('upper_dp', 4.2);
	run;

    data work.or2;
      length col5-col7 $200;
      set work.or1;

	  %*Set Subvar;
	  if substr(reverse(effect),1,1) = "N" then &subvar. = "N";
	  else if substr(reverse(effect),1,1) = "Y" then &subvar. = "Y";

	  %*Set Treatment and delete rows for ZS 5g v ZS 10g;
	  if substr(effect,1,14) = "TRT02PN 1 vs 2" then delete;
	  else if substr(effect,1,14) = "TRT02PN 1 vs 3" then trt02pn = 1;
	  else if substr(effect,1,14) = "TRT02PN 2 vs 3 " then trt02pn = 2;
  
      if not missing(oddsratioest) and oddsratioest <0.001 then col5 = "<0.001";
      else if not missing(oddsratioest) then col5 = put(oddsratioest, 5.2);

      if nmiss(lowercl, uppercl)=0 then col6 = cat("(", put(lowercl, &lower_dp.), ", ", put(uppercl, &upper_dp.), ")");
      else if missing(lowercl) and not missing(uppercl) then col6 = cat("(   NC, ", put(uppercl, &upper_dp.), ")");
      else if missing(lowercl) and not missing(uppercl) then col6 = cat("(", put(lowercl, &lower_dp.), ", NC)");
	  else if nmiss(lowercl, uppercl)=2 then col6 ="NC";

      if pvalue gt 0.999 then col7 = ">0.999";
      else if pvalue lt 0.001 then col7 = "<0.001";
      else col7 = put(pvalue,6.3);
    run;

    proc sort data = work.or2 out = work.or2s;
      by trt02pn &subvar.;
    run;

    %*Treatment-by-subgroup interaction: p-value only;
    data work.interaction1;
	  length col8 $200;
      set work.est1 (where =(variable=upcase("TRT02PN*&subvar") and classval0 in ("1" "2") and not missing(classval1)));

	  %*Set Treatment;
      if not missing(classval0) then trt02pn = input(classval0,2.);

	  %*Set Subvar;
	  if not missing(classval1) then &subvar = strip(classval1);

	  %*Format p-value;
	  if not missing(probchisq) and probchisq gt 0.999 then col8 = ">0.999";
      else if not missing(probchisq) and probchisq lt 0.001 then col8 = "<0.001";
      else if not missing(probchisq) then col8 = put(probchisq,6.3);
	  else if missing(probchisq) then col8 = "NC";
	
      keep trt02pn &subvar. col8;
    run;

    proc sort data = work.interaction1 out = work.interaction1s;
      by trt02pn &subvar.;
    run;

    %*Merge counts and odds ratio together;
	data work.stats1;
      merge work.n_counts2 
            work.or2s (keep = &subvar. trt02pn col:) 
            work.interaction1s;
      by trt02pn &subvar.;

	  if col8 = "NC" then do;
	    col5 = "  NC";
		col6 = "NC";
		col7 = "NC";
	  end;
    run;

    proc sort data = work.stats1 
               out = work.qc_&progname. (keep = col:);
      by ord1 trt02pn;
    run;

  %end;

  %*Compare to prod side*;
  proc compare base = output.&progname. compare = work.qc_&progname. errors warnings listall;
  run;

  %validation(program = &progname., level = output);

%mend qc_m_eff207;

