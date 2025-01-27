****************************************************************************************************
*
* Program Name      : t_tdo.sas
* Protocol / Study  : ORCHARD MLD ISE
* Type              : TABLE Macro
* Description       : To create outputs eff266g eff266h eff267a eff267c eff267d eff267e eff267f 
*                     eff267g eff267h eff267i eff267j eff267k eff267l eff269a eff269aa eff269ab 
*                     eff269ac eff269ad eff269c eff269d eff269e eff269f eff269g eff269h eff269i 
*                     eff269j eff269k eff269l eff269m eff269n eff269o eff269p eff269q eff269r 
*                     eff269s eff269sc eff269t eff269tc eff269u eff269uc eff269x eff269y eff269z 
*                     eff270b eff270c eff270d eff271a eff271b eff271c eff272ax eff272bx eff272cx 
*                     eff273ax eff273bx eff273cx
*
* Author            : Jess Sadler
* Date Completed    : 19JAN2022
* Input Dataset     : ADTTE
* Macros Used       : %RSTART %RSTOP 
* Reference Program Used : t_tdo.sas 
*
****************************************************************************************************
* Change History
*
* Changed by     : Jess Sadler 
* Date of Change : 27JAN2022 
* Reason         : Updated to add in suffix for proportion event free rows
*
* Changed by     : Jess Sadler 
* Date of Change : 09MAR2022 
* Reason         : Where clause added to subset tables further
*
* Changed by     : Jess Sadler 
* Date of Change : 13APR2022 
* Reason         : Updated after lengths updated across the studies in the ADaMs
*
* Changed by     : Jess Sadler 
* Date of Change : 02AUG2022 
* Reason         : Update row labels for total follow up to specify from birth
*                  Updated formatting to fit log rank test onto first page
*
* Changed by     : Jess Sadler 
* Date of Change : 18AUG2022 
* Reason         : Updated to add an option to turn off the Cox PH model in cases where assumptions do not hold
*
* Changed by     : Emily Foreman
* Date of Change : 24FEB2023 
* Reason         : Updated to remove footnote reference in last column for outputs not split by disease variant
*                  Updated summary stats presentation
*                  Added code for outputs created using shell SV9
*
* Changed by     : Emily Foreman
* Date of Change : 13MAR2023 
* Reason         : Updated for S3QC comments
*
* Changed by     : Prentis MacLean
* Date of Change : 04MAR2024
* Reason         : Removed P drive references

* Changed by     : Ryan Smith
* Date of Change : 04MAR2024
* Reason         : Updated to change all references to _name_ to be case dependent
********************************************************************************************************/

*======================================================================================================;
* PARAMETERS:
* DSIN = Dataset to be read in (defaults to ADTTE)
* SL = dataset used to get N numbers (defaults to ADSL)
* PARAMCD = Parameter of interest for treated and NHx subjects
* CRYOSIB = Parameter of interest for Cryo siblings 
* DISVAR = Disease variant for selecting the relevant Cryo siblings
* POP = Population flag (e.g. PSLIASFL)
* FN = Set to N if output does not require NH footnote (defaults to missing)
* WHERE = additional where clause if needed to subset subjects further
* ENDPOINT = endpoint of interest, defaults to sMFS and can be set to dosMFS, doOS, OS or OSalt
* OUT = output number used to create output dataset and read in for proc report 
* COX = defaults to N where the cox PH model is not included in the table if the assumptions do not hold
*       set to Y if the model is to be included in the table
*======================================================================================================;

%macro t_tdo(dsin     = adam.adtte,
             sl       = adam.adsl,
             paramcd  = ,
             cryosib  = ,
             disvar   = ,
             pop      = ,
       fn       = ,
             where    = ,
             endpoint = sMFS,
             out      = ,
             cox      = N);

* Clear datasets from work library;
proc datasets lib=work nolist memtype=data kill;
run; quit;

* Clear log before running;
dm 'log;clear;';run;quit;

%let trtopt = trta;
%let trtopt01 = trt01a;
%let trtoptn = trtan;
%let trtopt01n = trt01an;

%global ntrt1 ntrt2 ;

*** read in required data ***;  

data work.temp0;
  set &dsin. (where = (paramcd="&paramcd." and &pop="Y" &where.)); * JSS 9Mar22 Where clause added;
run;

* If there is cryo sibling data then include that;
%if &cryosib. ne  %then %do;
  data work.temp0b;
    set &dsin. (where = (paramcd="&cryosib." %if &disvar. ne %then %do; and csircdgp="&disvar." %end;));
    * Set the treatment group to Natural History as these are untreated siblings;
  * JSS 13Apr22 Removed format statement;
    &trtopt.="Natural History";
    &trtoptn.=1;
  run;

  data work.temp1;
    set work.temp0 
        work.temp0b;
  run;
%end;

%if &cryosib. = %then %do;

  data work.temp1;
    set work.temp0;
  run;

%end;

proc sort data = work.temp1;
  by &trtoptn.;
run;

*** calculate Big N numbers for column headers, this excludes Cryo siblings ***;
data work.adsl;
  set &sl. (where = (&pop.="Y"));
run;

proc freq data = work.adsl noprint ;
  tables &trtopt01n.*&trtopt01. / nocol norow nopercent out = work.count1 (rename =(&trtopt01n.=&trtoptn. count=N));
run;

*** create macro variables of N numbers in each subtype for use in column headers ***;
data _null_;
  set work.count1;
  call symputx("ntrt"||compress(put(&trtoptn.,1.)),N);
run;

%put &ntrt1.;
%put &ntrt2.;

*------------------------------------------------------------------------------------------------------------*;
* Rows for number of subjects included in the population and analysis, number of subject years of follow-up, *;
* number of subjects with an event                                       *;
*------------------------------------------------------------------------------------------------------------*; 

** get totals number of subjects for percentages **;
proc sql;
  create table work.totals as
  select count(distinct usubjid) as total, &trtoptn.
  from work.temp1
  group by &trtoptn.;
quit;

* These also form the 3rd line of the table for the version with Cryo siblings;
proc transpose data = work.totals out = work.row0_3 prefix=trt;
  var total;
  id &trtoptn.;
run;

** If there are cryo siblings then include the following rows for number of subjects in populations, number of 
additional subjects and total number of subjects in the analysis;
%if &cryosib. ne %then %do;

  ** the number of subjects in the population (excludes cryo siblings) **;
  proc transpose data = work.count1 out = work.row0_1 prefix=trt;
    var N ;
    id &trtoptn.;
  run;

  ** get the number of cryo siblings (additional subjects);
  proc freq data = work.temp0b noprint;
    tables &trtoptn.*&trtopt. / nocol norow nopercent out = work.count1b (rename =(count=ncryo));
  run; 

  proc transpose data = work.count1b out= work.row0_2 prefix=trt;
    var ncryo;
    id &trtoptn.;
  run;

  ** Set these together with the total number of subjects in the analysis and format data;
  data work.row0 (keep = col: section);
  length col1 col2 col3 col4 $200;
    set work.row0_1 
        work.row0_2 
        work.row0_3;

  if upcase(_name_)="N" then do;
      section = 1;
      col1="Number of subjects in population";
  end;
  else if upcase(_name_)="NCRYO" then do;
    section = 2;
      %if &fn.=N %then %do; col1="Number of additional subjects [1]"; %end;
    %else %do; col1 = "Number of additional subjects [2]"; %end;
  end;
  else if upcase(_name_)="TOTAL" then do;
      section = 3;
      col1="Total number of subjects in analysis";
  end;

    col2 = " ";

    if trt2 = . then col3 = "    0";
    else col3 = put(trt2,5.);

    if trt1 = . then col4 = "    0";
    else col4 = put(trt1,5.);
  run;

%end; 

* Number of subject-years of follow-up;
proc sql noprint;
  create table work.subyrs as 
  select distinct sum(aval) as subyrs, &trtoptn.
  from work.temp1
  group by &trtoptn.;
quit;

proc transpose data = work.subyrs out = work.row1_ prefix=trt;
  var subyrs;
  id &trtoptn.;
run;

** Format data;
* 2Aug22 Updated row labels to specify from birth;
data work.row1 (keep = col: section);
  length col1 col2 col3 col4 $200;
  set work.row1_;
  %if &endpoint.=doOS or &endpoint.=dosMFS %then %do;
    %if &cryosib. ne %then %do;
      %if &fn.=N %then %do; col1="Number of subject-years of follow-up from onset [2]"; %end;
      %else %do; col1="Number of subject-years of follow-up from onset [3]"; %end;
    %end;
    %if &cryosib. = %then %do;
      %if &fn.=N %then %do; col1="Number of subject-years of follow-up from onset [1]"; %end;
      %else %do; col1="Number of subject-years of follow-up from onset [2]"; %end;
    %end;
  %end;
  %else %do;
    %if &cryosib. ne %then %do;
      %if &fn.=N %then %do; col1="Number of subject-years of follow-up from birth [2]"; %end;
      %else %do; col1="Number of subject-years of follow-up from birth [3]"; %end;
    %end;
    %if &cryosib. = %then %do;
      %if &fn.=N %then %do; col1="Number of subject-years of follow-up from birth [1]"; %end;
      %else %do; col1="Number of subject-years of follow-up from birth [2]"; %end;
    %end;
  %end;
  col2 = " ";
  col3 = put(trt2,7.1);
  col4 = put(trt1,7.1);
  section = 4;
run;

** Number of events per treatment group ***;
proc sql noprint;
  create table work.events as 
  select count(distinct usubjid) as events, &trtoptn.
  from work.temp1 
  where cnsr = 0 
  group by &trtoptn.;
quit;

%let trtoptn=trtan;
*** Percentages of events and transpose ***;
data work.pevents (keep = events pevt &trtoptn.);
  length pevt $200.;
  merge work.events 
        work.totals;
  by &trtoptn.;
  *** number and percentage of events ***;
  if events = . then pevt = "    0";
  else if events ne . then pevt = cat(put(events, 5.), " (", strip(put(100*(events/total), 3.0)), "%)");
run;

proc transpose data = work.pevents out = work.row2_ prefix=trt;
  var pevt;
  id &trtoptn.;
run;

proc transpose data = work.pevents out = work.nevents (drop = _name_) prefix=trt;
  var events;
  id &trtoptn.;
run;

** Format data;
data work.row2;
  length col1 col2 col3 col4 $200;
  set work.row2_;
  %if &cryosib. ne %then %do;
    %if &fn.=N %then %do; 
      %if &endpoint.=OS %then %do; col1="Number (%) of subjects who died [3]"; %end;
    %else %do; col1="Number (%) of subjects with an event [3]"; %end;
  %end;
    %else %do; 
      %if &endpoint.=OS %then %do; col1="Number (%) of subjects who died [4]"; %end;
    %else %do; col1="Number (%) of subjects with an event [4]"; %end;
    %end;
  %end;
  %if &cryosib. = %then %do;
    %if &endpoint.=OS or &endpoint.=doOS %then %do; col1="Number (%) of subjects who died"; %end;
  %else %do; col1="Number (%) of subjects with an event"; %end;
  %end;
  col2 = " ";
  col3 = trt2;
  col4 = trt1;
  section = 5;
  keep col: section;
run;

*------------------------------------------------------------------------------------------------------------;
* Fishers exact test of proportion vs. natural history                                                       ;
*------------------------------------------------------------------------------------------------------------;

* Run the text if there are no cells with zero counts;
proc sql noprint;
  select sum(cnsr)
  into: cens1 trimmed 
  from work.temp1
  where cnsr = 1;

  select sum(cnsr)
  into: cens0 trimmed
  from work.temp1
  where cnsr = 0;
quit;

%put &=cens1 &=cens0;

%if %sysevalf(&cens0. ne %str()) and %sysevalf(&cens1. ne %str()) %then %do;

  ods select none;
  proc freq data = work.temp1;
    tables cnsr*trta / exact;
    ods output fishersexact = work.fishers; 
  run;
  ods select all;

  *Fishers exact p-value;
  data work.row3 (keep = col: section);
    length col1 col2 col3 $ 200;
    set work.fishers (where = (name1='XP2_FISH'));
    section=6; 
    col1="Fisher's exact test of proportion vs. Natural History";
    col2="p-value";
    if nvalue1 = . then col3 = "     NE";
    else if . < nvalue1 < 0.001 then col3="   <0.001";
    else if nvalue1 > 0.999 then col3="   >0.999";
    else if nvalue1 ne . then col3 = put(nvalue1, 9.3);
  run;

%end;

%else %do;

  data work.row3;
    length col1 col2 col3 $ 200;
    col1 = "Fisher's exact test of proportion vs. Natural History";
    col2 = "p-value";
    section = 6;
    col3 = "     NE";
  run;

%end;

*------------------------------------------------------------------------------------------------------------*;
* Median time to event, 95% CI, quartiles and min and max                            *;
*------------------------------------------------------------------------------------------------------------*; 

**get median 95% CI and quartiles of age at deaths**;
* PDF of log log and survival plots to check assumptions of proportional hazards;
ods _all_ close;
*PM 4MAR2024- changed P: to &__root.;
ods pdf file="&__root.\&client.\&project.\&study.\&re.\output\sandbox\&out._prophazard.PDF";
ods graphics on;
title "Proportional Hazards for output &out.";
footnote "&sysdate. &systime."; 
proc lifetest data = work.temp1 alpha = 0.05 plots=all;
  time aval*cnsr(1);
  strata &trtoptn. / test=none;
  ods output Quartiles = work.stats;
run;

* JSS 18Aug22 Added in case Cox model is not run;
%if &cox.=N %then %do;

  ods graphics off;
  ods pdf close;
  ods select all ;

%end;

data work.stats1 (keep = &trtoptn. median ci q1q3);
  length median ci q1q3 $50;
  set work.stats;

  if percent = 75 then do;
    if estimate = . then q3 = .;
  else q3 = estimate;
  end;

  if percent = 50 then do;
    if estimate = . then median = "     NE";
  else median = put(estimate, 7.1);

  if lowerlimit = . and upperlimit = . then ci = "    (NE, NE)";
  else if lowerlimit ne . and upperlimit = . then do;
      if lowerlimit lt 10 then ci = cat("   (", put(lowerlimit,3.1), ", NE)");
    else if lowerlimit ge 10 then ci = cat("  (", put(lowerlimit,4.1), ", NE)");
  end;
  else if lowerlimit = . and upperlimit ne . then ci = cat("   (NE, ", strip(put(upperlimit,4.1)), ")");
  else if lowerlimit ne . and upperlimit ne . then do;
      if lowerlimit lt 10 then ci = cat("   (", put(lowerlimit, 3.1), ", ", strip(put(upperlimit, 4.1)), ")");
    else if lowerlimit ge 10 then ci = cat("  (", put(lowerlimit, 4.1), ", ", strip(put(upperlimit, 4.1)), ")");
  end;
  end;

  if percent = 25 then do;
    if estimate = . and q3 = . then q1q3 = "     NE, NE";
  else if estimate ne . and q3 = . then q1q3 = cat(put(estimate, 7.1), ", NE");
  else if estimate = . and q3 ne . then q1q3 = cat("    NE, ", strip(put(q3, 6.1)));
  else if estimate ne . and q3 ne . then q1q3 = cat(put(estimate, 7.1), ", ", strip(put(q3, 6.1)));
  output;
  end;

  retain median ci q3 q1q3;
run;

proc transpose data = work.stats1 out = work.stats1_t prefix = trt;
  id &trtoptn.;
  var median ci q1q3;
run;

data work.stats2 (keep = ord col2 trt:);
  set work.stats1_t;
  if upcase(_name_) = "MEDIAN" then do;
    col2 = "Median";
  ord = 1;
  end;
  else if upcase(_name_) = "CI" then do;
    col2 = "95% CI";
  ord = 2;
  end;
  else if upcase(_name_) = "Q1Q3" then do;
    col2 = "Q1, Q3";
  ord = 3;
  end;
run;

* Median, min and max;
proc means data = work.temp1 (where=(cnsr=0)) MIN MAX MEDIAN N noprint;
  by &trtoptn.;
  var aval;
  output out = work.stats3 MIN=MIN MAX=MAX MEDIAN=MEDIAN N=N; 
run;

** Percentages of events and transpose ***;
data work.stats4;
  length mm medevt nevt $200.;
  set work.stats3;
  by &trtoptn.;

  ** min & max time to event ***;
  if min = . and max = . then mm="  NE";
  else if min ne . and max ne . then mm =put(min,7.1)||", "||strip(put(max,6.1));

  ** median time to event ***;
  if median = . then medevt="  NE";
  else if median ne . then medevt =put(median,7.1);

  ** n ***;
  if n = . then nevt="    0";
  else if n ne . then nevt=put(n,5.0);
run;

proc transpose data = work.stats4 out = work.stats4_t prefix=trt;
  var mm medevt nevt; 
  id &trtoptn.; 
run;

* Check if there are events;
data _null_;  
  dsid=open('stats4_t');
  chk1=varnum(dsid,'trt1');  
  chk2=varnum(dsid,'trt2');
  call symputx('chk1',chk1);  
  call symputx('chk2',chk2);
run;

%put &chk1. &chk2.; 

** sort out data for table **;
data work.row4a (keep=col: section ord);
  length col1 col2 col3 col4 $200.;
  set work.stats4_t;
  if upcase(_name_)="MM" then do; 
    col2 = "Min, Max"; 
    ord = 3;
  end;
  if upcase(_name_)="MEDEVT" then do; 
    col2 = "Median"; 
    ord = 2;
  end;
  if upcase(_name_)="NEVT" then do; 
    col2 = "n"; 
    ord = 1;
  end;

  %if &endpoint.=dosMFS %then %do;
    %if &cryosib. ne %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual time to event (years) for subjects with an event [4]"; %end;
        %else %do; col1="Actual time to event (years) for subjects with an event [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual time to event (years) for subjects with an event [2]"; %end;
        %else %do; col1="Actual time to event (years) for subjects with an event [3]"; %end;
    end;
    %end;
  %end;
  %else %if &endpoint.=doOS %then %do;
    %if &cryosib. ne %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual time to death (years) for subjects who died [4]"; %end;
        %else %do; col1="Actual time to death (years) for subjects who died [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual time to death (years) for subjects who died [2]"; %end;
        %else %do; col1="Actual time to death (years) for subjects who died [3]"; %end;
    end;
    %end;
  %end;
  %else %if &endpoint.=OS %then %do;
    %if &cryosib. ne %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual age at death (years) for subjects who died [4]"; %end;
        %else %do; col1="Actual age at death (years) for subjects who died [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual age at death (years) for subjects who died [2]"; %end;
        %else %do; col1="Actual age at death (years) for subjects who died [3]"; %end;
    end;
    %end;
  %end;
  %else %do;
    %if &cryosib. ne %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual age at event (years) for subjects with an event [4]"; %end;
        %else %do; col1="Actual age at event (years) for subjects with an event [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="n" then do;
        %if &fn.=N %then %do; col1="Actual age at event (years) for subjects with an event [2]"; %end;
        %else %do; col1="Actual age at event (years) for subjects with an event [3]"; %end;
    end;
    %end;
  %end;

  section=7.1;

  %if &chk2.>0 %then %do;
    col3 = trt2;
  %end;
  %if &chk1.>0 %then %do;
    col4 = trt1; 
  %end;
  %if &chk2.=0 %then %do;
    if upcase(_name_) in ("MM","MEDEVT") then col3="     NE"; 
    else if upcase(_name_) = "NEVT" then col3=put(0,5.0); 
  %end;
  %if &chk1.=0 %then %do;
    if upcase(_name_) in ("MM","MEDEVT") then col4="     NE"; 
    else if upcase(_name_) = "NEVT" then col4=put(0,5.0);  
  %end;
run;

** sort out data for table **;
data work.row4b (keep=col: section ord);
  length col1 col2 col3 col4 $200.;
  set work.stats2;
  %if &endpoint.=dosMFS %then %do;
    %if &cryosib. ne %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for time to event (years) [4]"; %end;
        %else %do; col1="Kaplan-Meier estimates for time to event (years) [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for time to event (years) [2]"; %end;
        %else %do; col1="Kaplan-Meier estimates for time to event (years) [3]"; %end;
    end;
    %end;
  %end;
  %else %if &endpoint.=doOS %then %do;
    %if &cryosib. ne %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for time to death (years) [4]"; %end;
        %else %do; col1="Kaplan-Meier estimates for time to death (years) [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for time to death (years) [2]"; %end;
        %else %do; col1="Kaplan-Meier estimates for time to death (years) [3]"; %end;
    end;
    %end;
  %end;
  %else %if &endpoint.=OS %then %do;
    %if &cryosib. ne %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for age at death (years) [4]"; %end;
        %else %do; col1="Kaplan-Meier estimates for age at death (years) [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for age at death (years) [2]"; %end;
        %else %do; col1="Kaplan-Meier estimates for age at death (years) [3]"; %end;
    end;
    %end;
  %end;
  %else %do;
    %if &cryosib. ne %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for age at event (years) [4]"; %end;
        %else %do; col1="Kaplan-Meier estimates for age at event (years) [5]"; %end;
    end;
    %end;  
    %if &cryosib. = %then %do;
      if col2="Median" then do;
        %if &fn.=N %then %do; col1="Kaplan-Meier estimates for age at event (years) [2]"; %end;
        %else %do; col1="Kaplan-Meier estimates for age at event (years) [3]"; %end;
    end;
    %end;
  %end;
  section=7.2;   
  col3 = trt2;
  col4 = trt1; 
run;

*------------------------------------------------------------------------------------------------------------*;
* Log rank test                                                                        *;
*------------------------------------------------------------------------------------------------------------*; 

** get estimate of unstratified log rank test**;
title "Unstratified log rank test for output &out.";
proc lifetest data = work.temp1 alpha = 0.05 method=km;
  time aval*cnsr(1);
  test &trtoptn.;
  ods output LogUniChiSq = work.lgranks;
run;  

** Get p-value for log-rank test **;
data work.row5 (keep = col: section);
  length col1 col2 trt2 col3 $200;
  set work.lgranks;
  section = 8;
  col1 = "Unstratified log rank test vs. Natural History";
  col2 = "p-value";
  if probchisq ne . and probchisq >= 0.001 then trt2 = put(round(probchisq, 0.001), 9.3);
  else if probchisq ne . and probchisq < 0.001 then trt2 = "   <0.001";
  if probchisq = . then trt2 = "     NE";
  col3 = trt2; 
run;

*------------------------------------------------------------------------------------------------------------*;
* Risk reduction - Cox PH model                                                        *;
*------------------------------------------------------------------------------------------------------------*;

* JSS 18Aug22 Added to avoid running the Cox PH model in case assumptions do not hold;
 
%if &cox.=Y %then %do;

  title "Cox proportional hazard test for output &out.";
  proc phreg data = work.temp1;
    class &trtopt. (ref = 'Natural History') / param = ref;
    model aval*cnsr(1)=&trtopt. / rl=wald alpha = 0.05 ties = efron maxiter=50;
    ods output parameterestimates = work.hazard;
  run;

  ods graphics off;
  ods pdf close;
  ods select all ;

  *format risk reduction table;
  data work.hazard2;
    length riskred ci rrll rrul $ 200;
    set work.hazard;
    *get risk reduction;
    if not missing(hazardratio) then do; 
      if hazardratio=0 then riskred = "     NE";
      else if 100/hazardratio>=99999.5 then riskred = put(100/hazardratio, best.)||"%";
      else riskred=put(100/hazardratio,7.1)||"%"; 
    end;
    if missing(hazardratio) then do; 
      riskred = "   NE";
    end;

    *get pvalue;
    if probchisq ne . and probchisq >= 0.001 then pvalue = put(probchisq, 9.3);
    else if probchisq ne . and probchisq < 0.001 then pvalue = "   <0.001";
    if probchisq = . then pvalue = "     NE";

    *upper and lower ci;
    if hruppercl in (0,.) then rrll='  (NE,';
    else if  0<100/hruppercl<9.5 then rrll='   ('||strip(put(100/hruppercl,8.1))||'%,';
    else if 9.5<=100/hruppercl<99.5 then rrll='  ('||strip(put(100/hruppercl,8.1))||'%,';
    else if 99.5<=100/hruppercl<999.5 then rrll=' ('||strip(put(100/hruppercl,8.1))||'%,'; 
    else if 999.5<=100/hruppercl then rrll='('||strip(put(100/hruppercl,8.1))||'%, ';
      
    if hrlowercl in (0,.) then rrul=' NE)'; 
    else if hrlowercl>0 then rrul=' ' ||strip(put(100/hrlowercl ,8.1))||'%)';

    if hruppercl in (0,.) and hrlowercl in (0,.) then ci="   NE";
    else if hruppercl not in (0,.) or hrlowercl not in (0,.) then ci=trim(rrll)||trim(rrul);
  run;

  *transpose risk hazard to get the correct format for the table;
  proc transpose data = work.hazard2 out = work.hazard3 prefix=colx;
    var riskred pvalue ci;
  run;

  data work.row6 (keep=col: section ord);
    length col1 col2 col3 $ 200;
    set work.hazard3;
    section = 9;
    %if &cryosib. ne %then %do;
    if upcase(_name_)="RISKRED" then do;
        %if &fn.=N %then %do; col1="Risk reduction vs. Natural History [5]"; %end;
        %else %do; col1="Risk reduction vs. Natural History [6]"; %end;
    end;
    %end;
    %if &cryosib. = %then %do;
      if upcase(_name_)="RISKRED" then do;
        %if &fn.=N %then %do; col1="Risk reduction vs. Natural History [3]"; %end;
        %else %do; col1="Risk reduction vs. Natural History [4]"; %end;
    end;
    %end;
    if upcase(_name_)="RISKRED" then do;
      col2="Risk reduction"; ord=1;
    end;
    else if upcase(_name_)="CI" then do;
      col2="95% CI"; ord=2;
    end;
    else if upcase(_name_)="pvalue" then do;
      col2="p-value"; ord=3;
    end;
    col3=colx1;
  run;

%end;

*------------------------------------------------------------------------------------------------------------*;
* Estimated proportion event free                                                      *;
*------------------------------------------------------------------------------------------------------------*;

** get maximum year in data**;
proc sql noprint;
  select max(ceil(aval)) 
  into: maxyr1 
  from work.temp1;
quit; 

** get estimates of survivial to 4,5,6, etc years **;
ods select none;
proc lifetest data = work.temp1 conftype=loglog timelist=(1 to &maxyr1. by 1) outs = work.yrsurv reduceout;
  time aval*cnsr(1);
  strata &trtoptn. /test=none;
run;
ods select all;

%let trtoptn=trtan;
** sort data for table by transposing and finding percs **;
data work.yrsurv1 (keep = timelist &trtoptn. disp1 disp2);
  length disp2 $200;
  merge work.yrsurv;
  by &trtoptn.;

  if survival = . then disp1 = "     NE";
  else if survival = 1 then disp1 = "    100";
  else if survival = 0 then disp1 = "      0";
  else disp1=put(round(survival*100,.1),7.1);

  if sdf_lcl = 1 and sdf_ucl = 1 then disp2 = "   (100, 100)";
  else if sdf_lcl = . and sdf_ucl ne . then disp2 = cat("    (NE, ", strip(put(100*sdf_ucl, 6.1)), ")");
  else if sdf_lcl ne . then do;
    if sdf_ucl = . then do;
      if 100*sdf_lcl <10 then disp2 = cat("   (", put(100*sdf_lcl, 3.1), ", NE)");
    else if 100*sdf_lcl ge 10 then disp2 = cat("  (", put(100*sdf_lcl, 4.1), ", NE)");
  end;
  else if sdf_ucl ne . then do;
      if 100*sdf_lcl <10 then disp2 = cat("   (", put(100*sdf_lcl, 3.1), ", ", strip(put(100*sdf_ucl, 6.1)), ")");
    else if 100*sdf_lcl ge 10 then disp2 = cat("  (", put(100*sdf_lcl, 4.1), ", ", strip(put(100*sdf_ucl, 6.1)), ")");
  end;
  end;
  else if sdf_lcl = . and sdf_ucl = . then disp2 = "    (NE, NE)";
run;

proc sort data = work.yrsurv1;
  by timelist;
run;

proc transpose data = work.yrsurv1 out = work.yrsurv2 prefix=trt;
  by timelist;
  var disp1 disp2;
  id &trtoptn.;
run;

%if &endpoint.=dosMFS or &endpoint.=doOS %then %do;

  data work.yrsurv5;
    set work.yrsurv2;
  run;

%end;
%else %do;

  * Get the last visit where both are 100% to get the start year;
  data work.startyr (keep=startyr one);
    set work.yrsurv2 (where = (lowcase(_name_)="disp1" and strip(trt1)="100" and strip(trt2)="100"));
    rename timelist=startyr;
    one=1;
  run;
  
  proc sort data = work.startyr out = work.s_startyr;
    by startyr one;
  run;
  
  data work.startyrb;
    set work.s_startyr;
    by one;
    if last.one;
  run;
  
  data work.yrsurv3;
    set work.yrsurv2;
    one=1;
  run;
  
  data work.yrsurv4;
    merge work.yrsurv3 
          work.startyrb;
    by one;
  run;
  
  data work.yrsurv5;
    set work.yrsurv4 (where = (timelist>=startyr or startyr=.));
  run;

%end;

data work.row7 (keep=col: section ord timelist);
  length col1 col2 col3 col4 $200;
  set work.yrsurv5;
  * JSS 27Jan22 Updated;
  if upcase(_name_) = "DISP1" then do; 
    %if &endpoint.=dosMFS %then %do;
    %if &cryosib. ne %then %do;
        %if &fn.=N %then %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year after onset [4]"; 
      else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years after onset [4]"; 
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year after onset [5]"; 
      else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years after onset [5]"; 
        %end;
    %end;
    %if &cryosib. = %then %do;
        %if &fn.=N %then %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year after onset [2]"; 
      else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years after onset [2]"; 
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year after onset [3]"; 
      else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years after onset [3]"; 
        %end;
    %end;
    %end;
    %else %if &endpoint.=doOS %then %do;
    %if &cryosib. ne %then %do;
        %if &fn.=N %then %do; 
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year after onset [4]"; 
          else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years after onset [4]";
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year after onset [5]"; 
          else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years after onset [5]";
        %end;
    %end;
    %if &cryosib. = %then %do;
        %if &fn.=N %then %do;
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year after onset [2]"; 
          else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years after onset [2]";
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year after onset [3]"; 
          else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years after onset [3]";
        %end;
    %end;
    %end;
    %else %if &endpoint.=OS %then %do;
    %if &cryosib. ne %then %do;
        %if &fn.=N %then %do; 
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year of age [4]"; 
      else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years of age [4]";
        %end;
        %else %do;  
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year of age [5]"; 
      else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years of age [5]";
        %end;
    %end;
    %if &cryosib. = %then %do;
        %if &fn.=N %then %do;  
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year of age [2]"; 
      else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years of age [2]";
        %end;
        %else %do;  
          if timelist = 1 then col1 = "Estimated proportion alive up to 1 year of age [3]"; 
      else col1 = "Estimated proportion alive up to "||strip(put(timelist,2.))||" years of age [3]";
        %end;
    %end;
    %end;
    %else %do;
    %if &cryosib. ne %then %do;
        %if &fn.=N %then %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year of age [4]";
          else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years of age [4]"; 
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year of age [5]";
          else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years of age [5]"; 
        %end;
    %end;
    %if &cryosib. = %then %do;
        %if &fn.=N %then %do;
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year of age [2]";
          else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years of age [2]"; 
        %end;
        %else %do; 
          if timelist = 1 then col1 = "Estimated proportion event-free up to 1 year of age [3]";
          else col1 = "Estimated proportion event-free up to "||strip(put(timelist,2.))||" years of age [3]"; 
        %end;
    %end;
    %end;
    col2 = "%"; ord=1; 
  end; 
  if upcase(_name_) = "DISP2" then do;  
    col2 = "95% CI"; ord= 2; 
  end; 
  col3 = trt2;
  col4 = trt1;
  %if &endpoint.=dosMFS or &endpoint.=doOS %then %do; section=timelist+10; %end;
  %else %do; section=timelist-startyr+10; %end;
run;

*------------------------------------------------------------------------------------------------------------*;
* Set up output dataset and report                                                   *;
*------------------------------------------------------------------------------------------------------------*;

**** set all data together and set up pagination;
* JSS 2Aug22 Updated to include log rank test on first page;
data work.all (keep=col: ord section pgn id);
  set %if &cryosib. ne %then %do; 
        work.row0 
      %end; 
      work.row1 
      work.row2 
      work.row3 
      work.row4a 
      work.row4b 
      work.row5 
      %if &cox.=Y %then %do; 
        work.row6 
      %end; 
      work.row7;
  %if &cox.=Y %then %do;
    %if &cryosib. ne %then %do;
      if section <=8 then pgn = 1;
      else if 9<=section<=14 then pgn = 2;
      else if 15<=section<=21 then pgn = 3;
      else if 22<=section<=28 then pgn = 4;
      else if 29<=section<=35 then pgn = 5;
      else pgn = 6;
    %end;
    %if &cryosib. = %then %do;
      if section <=9 then pgn = 1;
      else if 10<=section<=16 then pgn = 2;
      else if 17<=section<=23 then pgn = 3;
      else if 24<=section<=30 then pgn = 4;
      else if 31<=section<=37 then pgn = 5;
      else pgn = 6;
    %end;
  %end;
  %if &cox.=N %then %do;
    %if &cryosib. ne %then %do;
      if section <=8 then pgn = 1;
      else if 9<=section<=16 then pgn = 2;
      else if 17<=section<=23 then pgn = 3;
      else if 24<=section<=30 then pgn = 4;
      else if 31<=section<=37 then pgn = 5;
      else pgn = 6;
    %end;
    %if &cryosib. = %then %do;
      if section <=11 then pgn = 1;
      else if 12<=section<=18 then pgn = 2;
      else if 19<=section<=25 then pgn = 3;
      else if 26<=section<=32 then pgn = 4;
      else if 33<=section<=39 then pgn = 5;
      else pgn = 6;
    %end;
  %end;

  id = 1; *Create a variable for merging;

run;

*Merge on the number of events in each treatment group. If an arm has 0 events then set the summary stats to missing;
data work.nevents1;
  set work.nevents;
  id = 1;
run;

data work.all_1 (drop = id);
  merge work.all
        work.nevents1;
  by id;

  if section in (7.1, 7.2) then do;
    if trt2 = . then col3 = "";
  if trt1 = . then col4 = "";
  end;
run;
  
proc sort data = work.all_1 out=output.&out.; 
  by section ord; 
run; 

ods noresults;
%rstart;

* JSS 19Aug22 Updated to include log rank test on first page;
proc report data = output.&out. nowd split = "~" headline headskip nocenter missing spanrows;
  column PGN SECTION COL1 ORD COL2 COL3 COL4;

  define pgn     / order   order=internal noprint;
  define section / order   order=internal noprint;
  define col1    / display order=internal          ""                                style(header)={just=l} style(column)={just=l cellwidth=50% asis=on};
  define ord     / order   order=internal noprint;                                   
  define col2    / display order=internal          ""                                style(header)={just=l} style(column)={just=l cellwidth=12% asis=on vjust=bottom};
  define col3    /                                 "OTL-200~(N=&ntrt2.)"             style(header)={just=l} style(column)={just=l cellwidth=18% asis=on vjust=bottom};

  %if &fn.=N %then %do;
    define col4    /                               "Natural History~(N=&ntrt1.)"     style(header)={just=l} style(column)={just=l cellwidth=18% asis=on vjust=bottom};
  %end;
  %else %do;
    define col4    /                               "Natural History [1]~(N=&ntrt1.)" style(header)={just=l} style(column)={just=l cellwidth=18% asis=on vjust=bottom};
  %end;
  
  break after pgn / page; 

  compute after section / style={font_size=0.1pt};
    line " ";
  endcomp; 

run;

%rstop; 

%mend t_tdo; 

**** END OF MACRO PROGRAM ***;
