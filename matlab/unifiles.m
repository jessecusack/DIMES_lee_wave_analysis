% program UNIFILES -- Unify EM-APEX data into single file
%
% Grab the latest profiles and combine into one file. For this version
% (written after Philippines Straits deployment) try to build in flexibility
% between cruises so that any improvements can be used to go back and
% reprocess previous datasets (including cblast, eddies, climode and phil_ex).
%
% Name convention: this file "unifiles.m" resides in the ~/emapex/matlab
% directory and should recognize how to run by the experiment directory it is
% called from. Individual experiment versions are labeled unifiles_fran.m,
% etc. Note that these individual versions should become obsolete in the near
% future (I hope).
%
% If possible, data is taken from emavel output (which can be re-run if
% geomagnetic corrections need to be made, for example).

%clear all
clear outfile cx uni gps vit vel ctd efp

cruise_id = 'eddies'; % Finestructure data for Ledwell tracer release 2005
cruise_id = 'phil_iop'; % IOP cruises in Philippines Straits 2007-8
cruise_id = 'climode'; % Toole deployment in Gulf Stream 2007
cruise_id = 'phil_ex'; % Exploratory cruise in Philippines Straits 2007
cruise_id = 'canyon'; % Monterey Canyon 2008
cruise_id = 'cblast'; % Inaugural EM-APEX deployment in hurricane Frances 2004
cruise_id = 'philip';
cruise_id = 'dimes';

% Alternate method: determine which experiment from the calling
% directory. Defaults to the last one above if the calling directory doesn't
% contain one of the options. Cruise and directory name can be slightly
% different, so the following two cell arrays define the correspondence.
%dir_ids = {'frances' 'eddies' 'climode' 'philstraits'};
dir_ids = {'frances' 'eddies' 'climode' 'philstraits' 'canyon' 'dimes'};
%cruise_ids = {'cblast' 'eddies' 'climode' 'phil_ex'};
%cruise_ids = {'cblast' 'eddies' 'climode' 'phil_iop'};
%cruise_ids = {'cblast' 'eddies' 'climode' 'dec-2007'};
%cruise_ids = {'cblast' 'eddies' 'climode' 'philip'};
cruise_ids = {'cblast' 'eddies' 'climode' 'philip' 'canyon' 'dimes'};

if exist('monitor_on') & monitor_on
  cruise_id = monitor_cruise;
  disp(['Called by emamon. Running unifiles on ' cruise_id '.']);
else
wd = pwd;
whichcr = zeros(size(cruise_ids));
for iic = 1:length(dir_ids),
  if ~isempty(findstr(wd,dir_ids{iic})),
    whichcr(iic) = 1;
  end
end
if ~isempty(find(whichcr)),
  cruise_id = cruise_ids{find(whichcr)};
  disp(['Called from ' wd ' so running ' cruise_id ' version of unifiles.']);
else
  disp(['Defaulting to ' cruise_id]);
end
end

% Switches used for all experiments
ctdonly = 0; % Don't include velocity data
uponly = 0; % Don't include down casts
do_save = 1;
verbose = 0;
do_vbsplots = 0; % diagnostics from absolute velocity calc
do_plots_each = 0; % mostly buoyancy diagnostics
do_plots_after = 0; % all T,S,U,V
do_wplots = 0;
do_navplots = 0; % best to use this OR vbsplots, not both
pause_each = 0;
pause_after = 0;
use_v6 = 0;

% Constant parameters
pthresh = 15; % up-profile surfaced threshold
pthresh2 = 130; % down-profile surfaced threshold
% note: increased pthresh2 from 120 because of a specific case -- float 4814
% profile 253 is a down-prof that starts at ~125m and previous down+up are
% missing. Not sure what happened here, but try and treat this as a
% slightly-extended surface interval (only ~45 min total).

switch cruise_id
  case 'cblast'
    %fltids={'1633a','1634a','1636a'};
    fltids={'1633','1634','1636'};
    %matdir = '/Users/girton/frances/data/proc/';
    matdir = '/Users/girton/frances/data/dec/';
    velmatdir = '/Users/girton/frances/data/vel/mat/';
    rawdir = '/Users/girton/frances/data/raw/'; % or use dec files here?
    outdir = '/Users/girton/frances/data/info/';
    vitdir = matdir;

    date_deploy = {'20040831T130000','20040831T130000','20040831T130000'};
    date_recover = {'20041007T150000','','20041006T150000'};

    lastfastprof = [218 218 216]; % last "Fast Profiling Mode" up profile

    %Ndeps = 84; % uponly
    Ndeps = 85;

    Nefdeps = 183;

    %Nprofs = 443; % Frances up profiles
    Nprofs = 769; % all Frances profiles (up and down)

    magvar0 = -11.19; % CBLAST value
    ar = [21 25 -73 -68]; % Just enough to show floats

  case 'eddies'
    fltids={'1632a','1633a', '1633b', '1635a', '1635b', '1636a', '1636b', '1636c'};
    matdir = '/Users/girton/eddies/mat/';
    velmatdir = '/Users/girton/eddies/velmat/';
    outdir = '/Users/girton/eddies/info/';

    % Dates recorded in format 30 for datestr (taken from JBG event log file
    % which uses data recorded on board at the time of the event or from bridge
    % log afterward; i.e., not precise but probably good to within a couple of minutes)
    date_deploy = {'20050730T190500','20050720T171400','20050802T235100',...
	  '20050723T145500','20050802T224100','20050722T154500',...
	  '20050801T213500','20050901T202000'};
    date_recover = {'20050906T150000','20050801T231000','20050902T134500',...
	  '20050801T170000','','20050731T150000','20050831T210000',...
	  '20050905T200000'};

    Ndeps = 196; % for all EDDIES (and 1632a only)

    %Nefdeps = 58; % 1632a only
    Nefdeps = 600; % for all EDDIES

    %Nprofs = 801; % total number of up profiles (EDDIES)
    %Nprofs = 777; % all EDDIES up profiles (except 1636c) minus aborted ones
    %Nprofs = 337; % 1632a only
    Nprofs = 1578; % all EDDIES profiles (up and down) minus aborted ones

    magvar0 = -13.85; % EDDIES value
    ar = [29 31.5 -71 -66]; % covers just eddy and surrounding ship track

  case 'climode'
    fltids={'1633c','1636c'};
    fltids={'1633c','1636d'};
    % Confirm deploy and recover times with John Toole. All but 1633
    % deployment were estimated from emailed descriptions, not specific times.
    date_deploy = {'20070305T122500','20070213T213000'};
    date_recover = {'20070312T203000','20070317T123000'};

    %matdir = '/Users/girton/climode/decmat/'; % for ema2mat files
    %vitdir = '/Users/girton/climode/info/'; % for gps and vit files (all profs)
    %velmatdir = '/Users/girton/climode/velmat/';
    matdir = '/data/emapex/climode/dec/mat/';
    vitdir = matdir;
    velmatdir = '/data/emapex/climode/vel/mat/';
    %outdir = vitdir;
    outdir = '/Users/girton/climode/info/';

    Ndeps = 338; % (278 from 1633 alone)
    Nefdeps = 134; % (126 from 1633 alone)
    Nprofs = 562; % all CLIMODE (up and down) = 94+468
    magvar0 = -13.85; % stand-in from eddies (but overwritten by geomag calc later)

    ar = [34 40 293 306]; % All emapex trajectories

  case 'phil_ex'
    fltids={'1633d','1636d','1636e'};
    fltids={'1633d','1636e','1636f'};
    % Confirm deploy and recover times with Joe Martin. Specific times given
    % for deployments. 1633d recovery estimated from GPS.
    date_deploy = {'20070609T090300','20070616T152100','20070624T223500'};
    date_recover = {'20070613T083300','20070619T143000','20070628T094100'};

    % Joe's stated info on deployment positions...
    % 1633d: 10 30.037'  N, 121 45.005'  E (Bottom depth = 1377 m)
    % 1636e:  9  0.024'  N, 123 59.957'  E (Bottom depth = 1548 m)
    % 1636f: 11 34.3490' N, 121 24.9650' E
    %
    % 1633d recovery was by a fishing boat--possibly in midst of first GPS interval
    % 1636e (Bohol Sea) recovered by ship. Joe's notes say 1400 UTC but GPS
    % continues until 1630. Looks like trajectory kinks at 1430 so that's
    % most likely recovery time.
    % 1636f (north Mindoro) recovered by ship. Joe's notes say ship position
    % was 11 41.4'N, 121 9.6'E at 10:11 UTC 28 June 2007 (~30 min after
    % recovery). GPS quality dropped between 9:31 and 9:45, so 9:41 recovery sounds good.
    

    matdir = '/data/emapex/phil_ex/dec/mat/';
    velmatdir = '/data/emapex/phil_ex/vel/mat/';
    vitdir = matdir;
    outdir = '/Users/girton/philstraits/explor_info/';

    Ndeps = 248; % for all Philippines Exploratory (232,248,247)
    Nefdeps = 132; % (121,127,132)
    Nprofs = 158; % all Phil_ex profiles (60+46+52)
    magvar0 = -0.404; % standin Phil_ex value until variable mag params implemented

    ar = [8 17 118.5 126.5]; % All emapex trajectories

%  case 'phil_iop'
  case 'dec-2007'
    fltids={'3305a','1633e','3304a','1636e','3304b'};
    fltids={'1633d','1636e','1636f','1633e','1633f','1636g','3304a','3304b','3304c','3305a'}; % Expl and IOP
    fltids={'1633e','1633f','1636g','3304a','3304b','3304c','3305a'}; % IOP
    
    date_deploy = {'20071201T051200','20071202T143000','20071205T065300',...
	  '20071220T085111','20071221T093600'};
    date_recover = {'20071220T223000','20071216T093000','20071206T232300',...
	  '20071222T130000','20071225T000000'};

    % listing both Explor. and IOP
    date_deploy = {'20070609T090300','20070616T152100','20070624T223500',...
	  '20071202T143000','','20071220T085111','20071205T065300','20071221T093600','','20071201T051200'};
    date_recover = {'20070613T083300','20070619T163000','20070628T093000',...
	  '20071216T093000','','20071222T130000','20071206T232300','20071225T000000','','20071220T223000'};

    date_deploy = {'20071202T143000','20080113T143700','20071220T085111',...
	  '20071205T065300','20071221T093600','20080111T122000','20071201T051200'};
    date_recover = {'20071216T090300','20080115T080000','',...
	  '20071206T232300','20071222T161300','','20071220T223000'};

    
    % Bridge log deployment positions:
    % 3305a: 12 30.633' N, 121 46.095' E @ 0512z 12/01/07 (N Mindoro Strait)
    % 1633e: 11 31.181' N, 121 48.673' E @ 1430z 12/02/07 (Panay Strait)
    % 3304a: 11 16.647' N, 121 55.465' E @ 0653z 12/05/07 (Panay Strait)
    % 1636g: 12 49.824' N, 120 36.805' E @ 0851z 12/20/07 (log says 1635) (N Mindoro)
    % 3304b: 12 08.482' N, 122 19.692' E @ 0936z 12/21/07 (Sibuyan Sea--noisy EF)
    %
    % From Kevin Bartlett:
    % 3304c: Float #3304 was deployed at Station #13 (13.00695 N, 124.64433
    %   E) at 12:20, January 11, 2008, UTC.
    %   [later message mentioned 124 25.02 E launch but I think this is a typo]
    % 1633f: We deployed 1633 at CTD station 33 (09 50.00 N 125 00.00 E) at
    %   14:37 January 13, 2008 UTC. 
    %
    % Notes:
    % 3305a recovery was by WWF boat out of Palawan
    % 1633e recovered by Melville in Sulu Sea
    %   0903z 12/16/07 10 32.743N. 121 21.283E. Drifter 1633 on deck EW
    % 3304a recovery was by a fishing boat
    % 1636g was not recovered (vanished near N Mindoro sill early 12/31)
    %   last GPS fix at 00:16:26 31 Dec 2007; probably stuck on bottom at ~0100Z
    % 3304b recovered by Melville after short drift. Log:
    %   1600z 12/22/07 12 11.332N. 122 13.986E. Arrive vicinity EM-APEX BW
    %   1613z 12/22/07 12 11.651N. 122 13.751E. EM-APEX on deck BW
    % 1633f recovered by Melville when stuck on surface ("16:00 local time")
    % 3304c was not recovered (long drift in Pacific)
    
    matdir = '/data/emapex/dec-2007/dec/mat/';
    velmatdir = '/data/emapex/dec-2007/vel/mat/';
    vitdir = matdir;
    outdir = '/Users/girton/philstraits/info/';

    Ndeps = 248; % for all Philippines Exploratory
    Nefdeps = 132;
    Nprofs = 158; % all Phil_ex profiles
    magvar0 = -0.404; % standin Phil_ex value until variable mag params implemented

    ar = [8 17 118.5 126.5]; % All emapex trajectories
  case 'philip'
    fltids={'3305a','1633e','3304a','1636e','3304b'};
    fltids={'1633e','1633f','1636g','3304a','3304b','3304c','3305a'}; % IOP
    fltids={'4088a'}; % 4088 only
    fltids={'3305c','4088a'}; % 2009 IOP (allprofs09)
    fltids={'1633d','1636e','1636f','1633e','1633f',...
	  '1636g','3304a','3304b','3304c','3305a'}; % Expl and IOP (for allprofs08)
    fltids={'1633d','1636e','1636f','1633e','1633f',...
	  '1636g','3304a','3304b','3304c','3305a','3305c','4088a'}; % All PHILEX
    
    date_deploy = {'20071201T051200','20071202T143000','20071205T065300',...
	  '20071220T085111','20071221T093600'};
    date_recover = {'20071220T223000','20071216T093000','20071206T232300',...
	  '20071222T130000','20071225T000000'};

    % just 2008 IOP
    date_deploy = {'20071202T143000','20080113T143700','20071220T085111',...
	  '20071205T065300','20071221T093600','20080111T122000','20071201T051200'};
    date_recover = {'20071216T090300','20080115T080000','',...
	  '20071206T232300','20071222T161300','','20071220T223000'};

    % 4088 only
    date_deploy = {'20090323T120000'};
    date_recover = {''};

    % full 2009 IOP (allprofs09)
    date_deploy = {'20090210T171800','20090323T120000'};
    date_recover = {'20090318T043745',''};

    % all PhilEx 2007+2008 (allprofs08)
    date_deploy = {'20070609T090300','20070616T152100','20070624T223500',...
	  '20071202T143000','20080113T143700','20071220T085111',...
	  '20071205T065300','20071221T093600','20080111T122000','20071201T051200'};
    date_recover = {'20070613T083300','20070619T143000','20070628T094100',...
	  '20071216T090300','20080115T080000','',...
	  '20071206T232300','20071222T161300','','20071220T223000'};

    % all PhilEx 2007-9
    date_deploy = {'20070609T090300','20070616T152100','20070624T223500',...
	  '20071202T143000','20080113T143700','20071220T085111',...
	  '20071205T065300','20071221T093600','20080111T122000','20071201T051200','20090210T171800','20090323T120000'};
    date_recover = {'20070613T083300','20070619T143000','20070628T094100',...
	  '20071216T090300','20080115T080000','',...
	  '20071206T232300','20071222T161300','','20071220T223000','20090318T043745',''};

    % Notes on all PhilEx deployments (some duplicated above):
    
    % 1633d: 2007 Exploratory cruise (Joe Martin), Panay Strait, between sill
    % and MP2 location
    % deployed 10 30.037'  N, 121 45.005'  E (Bottom depth = 1377 m)
    % recovered by a fishing boat after 4 days of profiling, drifting northward.
    % Probably picked up in midst of first GPS interval

    % 1636e: 2007 Exploratory cruise (Joe Martin) (SW Bohol Sea)
    % [Dunlap's list calls this '1636d1']
    % deployed  9  0.024'  N, 123 59.957'  E (Bottom depth = 1548 m); 16 June
    % recovered by ship after 3 days (19 June) of SW drift.
    % Joe's notes say 1400 UTC but GPS
    % continues until 1630. Looks like trajectory kinks at 1430 so that's
    % most likely recovery time.

    % 1636f: 2007 Exploratory cruise (Joe Martin) (Mindoro Strait "mixing bowl")
    % [Dunlap's list calls this '1636d2' and only contains 6 profiles]
    % deployed 11 34.3490' N, 121 24.9650' E; 24 June 2007
    % recovered by ship after 3.5 days of NW drift. Joe's notes say ship position
    % was 11 41.4'N, 121 9.6'E at 10:11 UTC 28 June 2007 (~30 min after
    % recovery). GPS quality dropped between 9:31 and 9:45, so 9:41 recovery sounds good.

    % 3305a: Joint US-Philippines cruise (Helen Phillips)
    %  N Mindoro Strait and mixing bowl
    % deployed (according to bridge log)
    %  at 12 30.633' N, 121 46.095' E @ 0512z 12/01/07
    % launched between MP1 and sill, then drifted south past sill into mixing bowl
    % stuck on surface after 6 days, at 0300 UTC 12/7/07, then drifted SW
    % across shallow reef area
    % into Sulu Sea.
    % Recovery was by WWF boat out of Palawan on 12/20/07 (2230 UTC)

    % 1633e: Joint US-Philippines cruise (Helen Phillips)
    % Mixing bowl, Panay Strait, Sulu Sea
    % deployed 11 31.181' N, 121 48.673' E @ 1430z 12/02/07
    % deployed north of Panay sill, drifted south past sill and MP2 into Sulu Sea
    % 1633e recovered by Melville in Sulu Sea
    %   0903z 12/16/07 10 32.743N. 121 21.283E. Drifter 1633 on deck EW

    % 3304a: Joint US-Philippines cruise (Helen Phillips), Panay Strait
    % deployed 11 16.647' N, 121 55.465' E @ 0653z 12/05/07 (Panay Sill)
    % 3304a recovery was by a fishing boat after only 1.5 days of SE drift

    % 1636g: Mooring deployment cruise (James Girton)
    % [Dunlap's list calls this '1636e']
    %  N Mindoro along with MP1 mooring
    % deployed 12 49.824' N, 120 36.805' E @ 0851z 12/20/07 (log says 1635)
    % 1636g profiled for 11 days but was not recovered (vanished just before
    % reaching N Mindoro sill early 12/31)
    %   last GPS fix at 00:16:26 31 Dec 2007; probably stuck on bottom at ~0100Z

    % 3304b: Mooring deployment cruise (James Girton), Sibuyan Sea
    % deployed 12 08.482' N, 122 19.692' E @ 0936z 12/21/07
    % noisy EF--probably no good velocities
    % 3304b recovered by Melville after 1 day NE drift. Log:
    %   1600z 12/22/07 12 11.332N. 122 13.986E. Arrive vicinity EM-APEX BW
    %   1613z 12/22/07 12 11.651N. 122 13.751E. EM-APEX on deck BW
    %
    
    % 3304c: Regional survey cruise (Kevin Bartlett)
    %  Open Pacific NE of San Bernardino Strait
    % deployed at Station #13 (13.00695 N, 124.64433 E) at 12:20, January 11, 2008, UTC.
    %   [later message mentioned 124 25.02 E launch but I think this is a typo]
    % 3304c was not recovered (9.5 month drift in Pacific before batteries
    % ran low)
    % last profile 9/26/08 (hprof=354)
    % but data collection actually failed on 4/30/08
    % (hprof=281) because of a firmware error in wrapping data more than
    % twice around the 8 MB flash memory. GPS positions found by rtgps are still good after
    % this time, but possibly only once per profile. Other data will be
    % difficult to recover.

    % 1633f: 2008 Regional survey cruise (Kevin Bartlett), NE Bohol Sea
    % deployed at CTD station 33 (09 50.00 N 125 00.00 E) at 14:37 January 13, 2008 UTC. 
    % only made 2 profiles (1 round trip) before becoming stuck on surface
    % 1633f recovered by Melville when stuck on surface ("16:00 local time")
    %
    
    % 3305c: 2009 IOP Process cruise (Joe Martin) [3305b in Monterey Canyon, Aug 2008]
    % Launched in the middle of the Mindoro mixing bowl and drifted SW for 5 weeks
    % deployed 1718 UTC on 2/10/09
    % Recovered at end of Arnold Gordon's regional cruise just before
    % entering shallow reef area SW of strait (needed some tide-synchronous
    % profiling to work off the reef before ship arrived) 
    % last profile at 0420 UTC 3/18/09
    
    % 4088a: Deployed by Pierre Flament near K1 critical latitude during
    % final 2009 transit to Taiwan
    % first profile ended 3/23/09, 1445 UTC
    % left in SCS and survived for 9 months; not recovered
    % EF battery died completely in late Aug and was starting to look pretty
    % bad by late July (or even early July)
    % last profile 12/16/09
    
    matdir = '/data/emapex/philip/dec/mat/';
    velmatdir = '/data/emapex/philip/vel/mat/';
    vitdir = '/Users/girton/philstraits/info/';
    outdir = vitdir;
    webdir = '/Library/WebServer/Documents/philex';

    if uponly
      Nprofs = 472; % (30+23+26+103+1+79+14+9+140+47)
      Ndeps = 385; % (231,248,247,227,49,297,215,235,385,228)
      Nefdeps = 600;
    else
      Nprofs = 949; % all PhilEx profiles (60+46+52+206+2+158+28+18+281+98)
      Ndeps = 923; % all PhilEx (232,248,247,227,49,376,215,235,706,228)
      Nefdeps = 689; % all PhilEx (121,127,132,212,87,304,195,211,600,568)
    end
    magvar0 = -0.404; % standin Phil_ex value until variable mag params implemented
    
    ar = [10.3 13 120.2 122.4]; % Mindoro/Tablas Straits
    ar = [10.3 16 118 122.4]; % Mindoro/Tablas plus SCS (for allprofs09)
    ar = [8 17 118.5 126.5]; % All emapex trajectories
    ar = [8 17.5 117.5 126.4]; % All PHILEX emapex 2007-9

  case 'canyon'
    fltids={'3305'};

    date_deploy = {'20080822T000000'}; % from Point Sur
    date_recover = {'20080825T000000'}; % from John Martin

    %matdir = '/Users/girton/climode/decmat/'; % for ema2mat files
    %vitdir = '/Users/girton/climode/info/'; % for gps and vit files (all profs)
    %velmatdir = '/Users/girton/climode/velmat/';
    matdir = '/data/emapex/canyon/dec/mat/';
    vitdir = matdir;
    velmatdir = '/data/emapex/canyon/vel/mat/';
    outdir = vitdir;

    Ndeps = 100; % not determined yet
    Nefdeps = 80; % 
    Nprofs = 30; % 
    magvar0 = -13.85; % stand-in from eddies (but overwritten by geomag calc later)

    ar = [35 38 -122 -121]; % not confirme
  
  case 'dimes'

    % DIMES float catalog:
    %  3767,4086,4087 deployed by Byron Kilbourne from Feb 2009 Revelle US1
    %    tracer injection cruise. One failed after 1 month, other two after 9 months.
    %  4089,4090,4812,4813,4814,4815 deployed by Brian Guest from Feb 2010 Thompson US2
    %    tracer survey cruise. 4089 died immediately, 4815 lasted 1 year, other
    %    4 lasted 18-19 months.
    %  4976,4977 were UK purchases (Naveira Garabato and Smeed) and deployed
    %    by Byron Kilbourne (and UK folks) from James Cook UK2 cruise in Nov
    %    2010. Both were profiling continuously and lasted about 3.5 months.
    %  4594,4595 included FLBB (ONR PhilEx purchase) and were deployed by
    %    Brian Guest from James Clark Ross UK2.5 cruise in April 2011. 4594
    %    lasted 6 months and 4595 lasted 17 months (as of Sep
    %    2012). [although last call e-mail from 4595 seems to have been Feb
    %    2012--what is previous statement based on? (emamon gets Nov 2012 as
    %    last call.)]
    %  4596,4597 also included FLBB and were deployed by Luc Rainville from
    %    James Cook UK3 cruise in Feb 2012. Both were ballasted light but
    %    have made upper ocean profiles for 7 months (as of Sep 2012).
    
    % Some deployment info:
    %
    % 4596 launched at 2/17/12, 0510 GMT, at 54° 00' S 49° 44.3' W (Rainville)
    % 4597 launched at 2/17/12  0710 GMT, at 54° 00' S 50° 20.7' W (Rainville)
    
    % All DIMES floats so far (2009 and 2010)
    %fltids={'3767a' '4086a' '4087a' '4089a' '4090a' '4812a' '4813a' '4814a' '4815a'  };
    %date_deploy = {'20090209T233800' '20090203T184200' '20090216T192900' ...
    %	  '20100221T233400' '20100217T214200' '20100227T192800' '20100211T234400' '20100207T182300' '20100203T225500'  };
    %date_recover = {'' '' '' '' '' '' '' '' ''};

    % Just 2009 floats
    %fltids={'3767a' '4086a' '4087a'};
    %date_deploy = {'20090209T233800' '20090203T184200' '20090216T192900'};
    %date_recover = {'' '' ''};

    % Just 2010 floats
    %fltids={'4089a' '4090a' '4812a' '4813a' '4814a' '4815a'};
    %date_deploy = {'20100221T233400' '20100217T214200' '20100227T192800' '20100211T234400' '20100207T182300' '20100203T225500'  };
    %date_recover = {'' '' '' '' '' ''};

    % Just 2011 floats
    %fltids={'4976a' '4977a' '4594a' '4595a'};
    %date_deploy = {'20101230T080000' '20101230T080000' '20110413T110000' '20110415T070000'};
    %date_recover = {'' '' '' ''};

    % 2010 and 2011 and 2012 floats
    %fltids={'4089a' '4090a' '4812a' '4813a' '4814a' '4815a' '4976a' '4977a' ...
%	  '4594a' '4595a' '4596a' '4597a'};
    %date_deploy = {'20100221T233400' '20100217T214200' '20100227T192800' ...
%	  '20100211T234400' '20100207T182300' '20100203T225500' ...
%	  '20101230T080000' '20101230T080000' ...
%	  '20110413T110000' '20110415T070000' ... 
%	  '20120217T051000' '20120217T071000'};
   %date_recover = {'' '' '' '' '' '' '' '' '' '' '' ''};
    
    % just UK floats
    %fltids={'4976a' '4977a'};
    %date_deploy = {'20101230T080000' '20101230T080000'};
    %date_recover = {'' ''};

    % just FLBBx floats (also current set as of fall 2012)
    %fltids={'4594a' '4595a' '4596a' '4597a'};
    %date_deploy = {'20110413T110000' '20110415T070000' '20120217T051000' '20120217T071000'};
    %date_recover = {'' '' '' ''};

    % floats deployed from L.M. Gould, Sept 2013
    %fltids={'6478a' '6480a' '6481a' '6625a' '6626a'};
    %date_deploy = {'20130916T080945'  '20130915T151720'  '20130915T201800'...
	%  '20130915T050335' '20130916T003830'};
    %lat_deploy = [-60-30.33/60 -58-00.68/60 -58-48.4/60 -56-30.0/60 -59-20.06/60];
    %lon_deploy = [-64-59.7/60 -65-00.14/60 -64-59.89/60 -64-59.6/60 -65-00.0/60];
    %date_recover = {'' '' '' '' ''};

    % all floats active in 2013 (including FLBBs and 5 deployed from L.M. Gould, Sept 2013)
    fltids={'4594a' '4595a' '4596a' '4597a' '6478a' '6480a' '6481a' '6625a' '6626a'};
    date_deploy = {'20110413T110000' '20110415T070000' '20120217T051000' ...
	  '20120217T071000' '20130916T080945'  '20130915T151720' ...
	  '20130915T201800' '20130915T050335' '20130916T003830'};
    lat_deploy = [NaN NaN -54 -54 -60-30.33/60 -58-00.68/60 -58-48.4/60 ...
	  -56-30.0/60 -59-20.06/60];
    lon_deploy = [NaN NaN -49-44.3/60 -50-20.7/60 -64-59.7/60 -65-00.14/60 ...
	  -64-59.89/60 -64-59.6/60 -65-00.0/60];
    date_recover = {'' '' '' '' '' '' '' '' ''};
    
    %matdir = '/Users/girton/climode/decmat/'; % for ema2mat files
    %vitdir = '/Users/girton/climode/info/'; % for gps and vit files (all profs)
    %velmatdir = '/Users/girton/climode/velmat/';
    %matdir = '/data/emapex/dimes/dec/mat/';
    matdir = '/emapex/dimes/dec/mat/';
    %velmatdir = '/data/emapex/dimes/vel/mat/';
    velmatdir = '/emapex/dimes/vel/mat/';
    vitdir = '/Users/girton/dimes/info/';
    outdir = vitdir;
    webdir = '/Library/WebServer/Documents/dimes';

    %Ndeps = 936; % for 2009: 936/923/920
    %Ndeps = 944; % for 2010: 794/943/920/916/930/944
    %Ndeps = 926; % for 2011: 765/762/926/687
    %Ndeps = 765; % for UK: 765/762
    %Ndeps = 926; % for FLBB: 926/687/453/577
    Ndeps = 999; % for Gould 2013

    %Nefdeps = 677; % for 2009: 625/625/677
    %Nefdeps = 692; % for 2010: 625/656/692/643/634/640
    %Nefdeps = 623; % for 2011: 482/480/623/582
    %Nefdeps = 482; % for UK: 482/480
    %Nefdeps = 683; % for FLBB: 623/582/332/683
    Nefdeps = 787; % for Gould 2013

    %Nprofs = 516; % UK floats as of Feb 2011
    %Nprofs = 712; % 2009 floats: 320+298+94
    %Nprofs = 1435; % 2010 floats: 4+260+291+300+258+322
    %Nprofs = 1789; % 2011 floats: 498+494+409+388
    %Nprofs = 992; % UK floats: 498+494
    %Nprofs = 1523; % FLBB floats: 409+388+439+287
    Nprofs = 100; % for Gould 2013

    magvar0 = NaN; % stand-in from eddies (but overwritten by geomag calc later)

    %ar = [-59 -58 -107.5 -105.5];
    %ar = [-64 -55 -107.5 -60];
    %ar = [-68 -40 -115 -25]; % largest field of view for full suite of floats
    %ar = [-65 -35 -90 -20]; % shift westward for 2011-2012 floats
    ar = [-65 -35 -90 -20]; % shift westward for 2011-2012 floats

end

if exist('monitor_on') & monitor_on
  Nprofs = nprofs_mon;
  disp(['Nprofs = ' num2str(nprofs_mon) ' before removing skipped profiles.']);
end

if ctdonly
  if use_v6
    outfile = [outdir 'allctd-v6.mat'];
  else
    outfile = [outdir 'allctd.mat'];
  end
else % if ctdonly
  if uponly
    outfile = [outdir 'allprofs_up.mat'];
  else
    outfile = [outdir 'allprofs.mat'];
  end
end

%if exist([outdir 'allprofs.mat']),
%uni = load([outdir 'allprofs.mat']);
if exist(outfile),
  uni = load(outfile);
  if ~isfield(uni,'Nskip'),
    uni.Nskip = 0;
  end
  Nprofs = Nprofs - uni.Nskip;
  if length(find(isfinite(uni.hpid)))==Nprofs,
    gotuni = 1;
    disp('using previous unifiles output');
  else
    gotuni = 0;
    disp(['loaded previous unifiles output but not long enough--next run' ...
	  ' will be better!']);
  end
else
  gotuni = 0;
end

nmoverlap = 5;
ar2 = ar + nmoverlap/60*[-1 1 -1 1];
[topo tlat tlon] = extract_1m(ar2); % new 1-min S+S bathy

% Initialize variables

[T,S,P,UTC] = deal(repmat(NaN,Ndeps,Nprofs));
[utc_gps,hpid,ssid,vbsid,lat_gps,lon_gps,flid,depl,dnup,ind,...
      utc_dep,utc_rec,utc_up,utc_down,maxp,botdep,...
      ubs,vbs,u_sfc,v_sfc,u_gps,v_gps,gps_dist,apf9_terr,nvel] =...
    deal(repmat(NaN,1,Nprofs));
[first_apf9time,median_gpstime,median_apf9time,...
      lat_up,lon_up,lat_down,lon_down,surfaced,argo_mode] =...
    deal(repmat(NaN,1,Nprofs));
[Wp_ca,ppos_ca,P_ca,mlt_ca] = deal(repmat(NaN,Ndeps-1,Nprofs));
got_gps = repmat(0,1,Nprofs);
Nskip = 0;

if ~ctdonly
[U,V,U1,V1,U2,V2,Pef,UTCef] = deal(repmat(NaN,Nefdeps,Nprofs));
%[E1sdev,E2sdev,V1woW,V2woW] = deal(repmat(NaN,Nefdeps,Nprofs));
[Verr1,Verr2,V1woW,V2woW,Vchan] = deal(repmat(NaN,Nefdeps,Nprofs));
%magvar = repmat(magvar0,1,Nprofs);
[magvar,fh,fz] = deal(repmat(NaN,1,Nprofs));
[Wp,Wr,Wf,ppos] = deal(repmat(NaN,Nefdeps,Nprofs));
end

uxt0 = datenum(1970,1,1,0,0,0);
ii = 1; % index from 1 to Nprofs

% Loop over all instruments
for fltidind=1:length(fltids),
%for fltidind=3,
%for fltidind=6:length(fltids),
  profnum = 1;
  vbsnum = 0; % first profile is number 0 since it won't have easy vbs
             % calculation data (until launch GPS is added in)
  %ind(ii) = ii;
  fltid = fltids{fltidind}; % float/deployment string

  switch fltid,
    % Correct for different labelling in John Dunlap's varinit.m database
    case '1636e',
      rfltid = '1636d1';
    case '1636f',
      rfltid = '1636d2';
    case '1636g',
      rfltid = '1636e';
    otherwise
      rfltid = fltid;
  end
  
  fid = str2num(fltid(1:4)); % float number
  if length(fltid)>4
    depchar = fltid(5); % deployment char
  else
    depchar = 'a';
  end
  switch depchar,
    case 'a',
      deploy = 1;
    case 'b',
      deploy = 2;
    case 'c',
      deploy = 3;
    case 'd',
      deploy = 4;
    case 'e',
      deploy = 5;
    case 'f',
      deploy = 6;
    case 'g',
      deploy = 7;
    otherwise
      error('depchar should be a, b, c, d, e, f, or g');
  end

  ddep_s = date_deploy{fltidind};
  ddep = datenum(ddep_s,'yyyymmddTHHMMSS');
  drec_s = date_recover{fltidind};
  if ~isempty(drec_s),
    drec = datenum(drec_s,'yyyymmddTHHMMSS');
  else
    drec = NaN;
  end
  
  % Loop over all CTD profiles (assuming that there are no EFP or other files
  % without a corresponding CTD)
  %
  % Finally (12/21/06) realized that velmat files (generated by emavel.m)
  % have essentially the same file structure as CBLAST ema files. So use
  % those unless a problem appears...
  
  if strcmp(cruise_id,'cblast'),
    %D = dir([matdir 'ema-' fltid(1:4) '-*.mat']);
    D = dir([velmatdir 'ema-' fltid '/*-vel-*.mat']);
    nD = size(D,1);
    gps = load([matdir 'ema-' fltid '/ema-' fltid '-gps.mat']);
  elseif strcmp(cruise_id,'dimes')|strcmp(cruise_id,'philip'),
    D = dir([velmatdir 'ema-' fltid '/*-vel.mat']);
    nD = size(D,1);
    %gps = load([matdir 'ema-' fltid '/ema-' fltid '-gps.mat']);
      vit = load([vitdir 'ema-' rfltid '-vit.mat']);
      gps = load([vitdir 'ema-' rfltid '-gps.mat']);
  else
    %D = dir([matdir 'ema-' fltid '/*-ctd-*.mat']);
    D = dir([velmatdir 'ema-' fltid '/*-vel-*.mat']);
    nD = size(D,1);
    if strcmp(cruise_id,'eddies')|strcmp(vitdir,matdir)
      vit = load([matdir 'ema-' fltid '/ema-' fltid '-vit.mat']);
      gps = load([matdir 'ema-' fltid '/ema-' fltid '-gps.mat']);
    else
      vit = load([vitdir 'ema-' rfltid '-vit.mat']);
      gps = load([vitdir 'ema-' rfltid '-gps.mat']);
    end
  end
  
  if verbose
    % Compare UTC of first profile with deployment time...
%    if strcmp(cruise_id,'cblast'),
%      emaname = ['ema-' fltid(1:4) '-0001.mat'];
      %vel = load([matdir emaname]);
%      vel = load([velmatdir 'ema-' fltid '/' emaname]);
%    else
    if strcmp(cruise_id,'dimes')|strcmp(cruise_id,'philip'),
      emaname = ['ema-' rfltid '-0001-vel.mat'];
    else
      emaname = ['ema-' rfltid '-vel-0001.mat'];
    end
    vel = load([velmatdir 'ema-' fltid '/' emaname]);
%   end
    disp([fltid ' launched ' datestr(ddep)]);
    for i = 1:3,
      %disp([datestr(ctd.UXT(i)/86400 + uxt0) ' p=' num2str(ctd.P(i))]);
      disp([datestr(vel.ctd_mlt(i)) ' p=' num2str(vel.Pctd(i))]);
    end
  end

  if ~ctdonly
  Uacc = []; Vacc = []; Timacc = []; Pacc = [];
  maxefvals = 0;
  end
  [utc_desc,lat_desc,lon_desc] = deal(NaN);
  [utc_desc_vbs,lat_desc_vbs,lon_desc_vbs] = deal(NaN);
  hpids = [];
  maxvals = 0;
  nmiss = 0; upinds = []; dninds = [];
  for iD = 1:nD,
%  for iD = 220,
%  for iD = 244, % misbehaving near-end profile of 1635b
    % read each CTD filename and save the hpid (if it is even, for uponly)
    %[aa,n] = sscanf(D(iD).name,'ema-%4d%1s-ctd-%4d');
      if strcmp(cruise_id,'cblast'),
	%[aa,n] = sscanf(D(iD).name,'ema-%4d-%4d');
	[aa,n] = sscanf(D(iD).name,'ema-%4d-vel-%4d');
	hprof = aa(2);
      elseif strcmp(cruise_id,'dimes')|strcmp(cruise_id,'philip'),
	if strcmp(fltid,'1636e')|strcmp(fltid,'1636f'),
	  [aa,n] = sscanf(D(iD).name,'ema-%4d%2s-%4d-vel');
	  hprof = aa(4);
	else
	  [aa,n] = sscanf(D(iD).name,'ema-%4d%1s-%4d-vel');
	  hprof = aa(3);
	end
      else
	[aa,n] = sscanf(D(iD).name,'ema-%4d%1s-vel-%4d');
	hprof = aa(3);
      end
    
    
    if ~uponly | mod(hprof,2)==0,
      hpids = [hpids hprof]; % hp listing for one deployment
      
      if strcmp(cruise_id,'dimes')|strcmp(cruise_id,'philip'),
	%emaname = ['ema-' fltid(1:4) '-' sprintf('%04d',hprof) '.mat'];
	velname = ['ema-' rfltid '-' sprintf('%04d',hprof) '-vel.mat'];
	%vel = load([matdir emaname]);
	vel = load([velmatdir 'ema-' fltid '/' velname]);
	%ctd = vel.CTD;
	% CTD file should not need to be loaded because all necessary info
	% has been transferred by emavel.m to velmat file!
      else
	%ctdname = ['ema-' fltid '-ctd-' sprintf('%04d',hprof) '.mat'];
	%ctd = load([matdir 'ema-' fltid '/' ctdname]);
	velname = ['ema-' rfltid '-vel-' sprintf('%04d',hprof) '.mat'];
	vel = load([velmatdir 'ema-' fltid '/' velname]);
      end
  
      % First time just loop through profiles and list number of hprofs per
      % deployment, also max # of pressure levels

      % Save profiles
      %npvals = length(ctd.P);
      npvals = length(vel.Pctd);
      if npvals<=2
%	if verbose
        disp(['Skipped float ' fltid ' half-profile# ' num2str(hprof) ' because only '...
	      num2str(npvals) ' data point.']);
	Nskip = Nskip+1;
%	end
      else
	
	% Now this is a real profile...
	if verbose
	  disp(['Starting to process profile ' num2str(hprof)]);
	end

	ind(ii) = ii;
	ssid(ii) = profnum;
	vbsid(ii) = vbsnum;
	%pmax(ii) = max(ctd.P);
	depl(ii) = deploy;
	flid(ii) = fid;
	hpid(ii) = hprof;
	dnup(ii) = 2-mod(hprof,2);
	
	if dnup(ii)==2 & min(vel.Pctd)<pthresh,
	  % up profile
	  surfaced(ii) = 1;
%	elseif dnup(ii)==1 & (hprof==1 | surfaced(ii-1)) & min(vel.Pctd)<pthresh2,
	elseif dnup(ii)==1 & (hprof==1 | surfaced(ii-1)),
	  % down profile
	  % Removed minimum pressure requirement because of float 3767
	  % profiles 123-125 which apparently reached the surface but don't
	  % have near-surface data. Is this ok?
	  surfaced(ii) = 1;
	else
	  surfaced(ii) = 0;
	end
	
	if dnup(ii)==2 & (ii==1 | dnup(ii-1)==2),
	  argo_mode(ii) = 1;
	else
	  argo_mode(ii) = 0;
	end
	
      %P(1:npvals,ii) = ctd.P;
      %T(1:npvals,ii) = ctd.T;
      %S(1:npvals,ii) = ctd.S;
      %UTC(1:npvals,ii) = ctd.UXT/86400 + uxt0;
      P(1:npvals,ii) = vel.Pctd;
      T(1:npvals,ii) = vel.T;
      S(1:npvals,ii) = vel.S;
      UTC(1:npvals,ii) = vel.ctd_mlt;

      magvar(ii) = vel.magvar;
      fh(ii) = vel.fh;
      fz(ii) = vel.fz;
      
      maxp(ii) = max(vel.Pctd);

      % UTC variable comes from ctd file, which probably means APF9 real time
      % clock that is set at start of deployment (at least that's my guess)
      
      if strcmp(cruise_id,'cblast'),
	if 1
	  % For CBLAST data, could apply a time correction here using 1.21
	  % sec/day drift diagnosed from comparison with GPS.
	  d0 = datenum(2004,8,1,12,0,0);
	  dcor = (vel.ctd_mlt(1)-d0)*1.21/86400;
	  vel.ctd_mlt = vel.ctd_mlt + dcor;
	end
      end
      
      if ~ctdonly
	if strcmp(cruise_id,'cblast'),
	if 1
	  % apply time correction to efp
	  dcor = (vel.efp_mlt(1)-d0)*1.21/86400;
	  vel.efp_mlt = vel.efp_mlt + dcor;
	end
      end

      if 0
      if strcmp(cruise_id,'cblast'),
	D2 = dir([matdir 'ema-' fltid '/' 'ema-' fltid '-efp-*-' sprintf('%04d',hprof) '.mat']);
	efpname = D2.name;
	else
	efpname = ['ema-' rfltid '-efp-' sprintf('%04d',hprof) '.mat'];
      end
      efp = load([matdir 'ema-' fltid '/' efpname]);
      end
      
	% Mask out EF grid pressure when CTD time gap is 30min or
	% greater. Also mask W computations.
	toolong = 30*60; % 1/2 hour
	lints = find(diff(vel.ctd_mlt)*86400>toolong);
	for li = lints,
	  mlt1 = vel.ctd_mlt(li);
	  mlt2 = vel.ctd_mlt(li+1);
	  bd = find(vel.efp_mlt>mlt1&vel.efp_mlt<mlt2);
	  vel.Pef(bd) = NaN*bd;
	  vel.Wef(bd) = NaN*bd;
	end
	
	% Note: these variables aren't saved in the efp file but only in the
	% emadec output. Need to reconstruct them using pressure
	% interpolation and velocity computation from another file...
	% or instead, use velmat file.
	gd = find(isfinite(vel.Pef)&isfinite(vel.efp_mlt)&isfinite(vel.Wef));
	bdw = find(isfinite(vel.Pef)&isfinite(vel.efp_mlt)&~isfinite(vel.Wef));
	if verbose &~isempty(bdw),
	  warning([num2str(length(bdw)) ' bad Wef velocity values not included.']);
	end
	nefpvals = length(vel.Pef(gd));
	Pef(1:nefpvals,ii) = vel.Pef(gd);
	UTCef(1:nefpvals,ii) = vel.efp_mlt(gd);

	% Select preferred velocity channel to save as U,V based on "verr" screening:
	% - first, remove all velocities with Verr>=0.1
	% - second, select quieter channel if Verr is half that on noiser channel
	% - finally, average the two if Verr has comparable levels
	% - set Vchan flag to 1 or 2 if only one channel used, 0 if both, and
	%   NaN if neither
	vecut = 0.1; % maximum single-point error
	vrcut = 0.5; % minimum error ratio for averaging
	
	% Correct a couple of mis-wired floats (until John builds it into his processing)
%	if flid(ii)==6625 | flid(ii)==6626,
%	  vel.V1 = -vel.V1;
%	  vel.U2 = -vel.U2;
%	end

	
	% Create temporary u,v vectors
	U1t = vel.U1(gd); U2t = vel.U2(gd); % row vectors
	V1t = vel.V1(gd); V2t = vel.V2(gd);

	% cut out values with high Verr
	bd1 = find(vel.Verr1(gd)>=vecut);
	U1t(bd1) = NaN; V1t(bd1) = NaN;
	bd2 = find(vel.Verr2(gd)>=vecut);
	U2t(bd2) = NaN; V2t(bd2) = NaN;
	
	% also remove known bad channels:
	% all CH2 for hpid>=48 on EMA 4813
	if flid(ii)==4813 & hpid(ii)>=48,
	  U2t = NaN*U2t; V2t = NaN*V2t;
	end
	
	warning off
	ve_ratio = vel.Verr1(gd)./vel.Verr2(gd);
	warning on
	
	ch1lo = find(ve_ratio < vrcut);
	U2t(ch1lo) = NaN; V2t(ch1lo) = NaN;
	ch2lo = find(ve_ratio > 1/vrcut);
	U1t(ch2lo) = NaN; V1t(ch2lo) = NaN;
	
	Ut = mymean([U1t; U2t]); Vt = mymean([V1t; V2t]);
	Vch = Ut*0;
	ch1 = find(isfinite(U1t)&~isfinite(U2t)); Vch(ch1) = 1;
	ch2 = find(~isfinite(U1t)&isfinite(U2t)); Vch(ch2) = 2;
	Vchan(1:nefpvals,ii) = Vch;
	
%	if flid(ii)==1633 & hpid(ii)==8,
%	  error('stopping to check top U,V value.')
%	end
	
	if verbose
	  % report how many velocity values started with, how many removed by
	  % Verr criterion, how many remaining and on which channel
	  disp(['From ' num2str(nefpvals) ' EF samples, dropped ' ...
		num2str(length(bd1)) ' from CH1 and ' num2str(length(bd2)) ...
		' from CH2 due to Verr cutoff.']);
	  disp(['Left with ' num2str(length(find(isfinite(Ut)))) ...
		' good samples, ' num2str(length(ch1)) ' CH1 only and ' ...
		num2str(length(ch2)) ' CH2 only.']);
	end
	
	
	[Ug,Vg] = vrot(Ut',Vt',magvar(ii));
	%[Ug,Vg] = vrot(vel.U1(gd)',vel.V1(gd)',magvar(ii));
	%[Ug,Vg] = vrot(ema.U2',ema.V2',magvar(ii));
	U(1:nefpvals,ii) = Ug;
	V(1:nefpvals,ii) = Vg;

	nvel(ii) = length(find(isfinite(Ug)));
	if verbose
	  disp(['nvel(ii) = ' num2str(nvel(ii))]);
	end

	U1(1:nefpvals,ii) = vel.U1(gd);
	V1(1:nefpvals,ii) = vel.V1(gd);
	U2(1:nefpvals,ii) = vel.U2(gd);
	V2(1:nefpvals,ii) = vel.V2(gd);

%	E1sdev(1:nefpvals,ii) = vel.e1sdev;
%	E2sdev(1:nefpvals,ii) = vel.e2sdev;
	Verr1(1:nefpvals,ii) = vel.Verr1(gd);
	Verr2(1:nefpvals,ii) = vel.Verr2(gd);
	V1woW(1:nefpvals,ii) = vel.V1woW(gd);
	V2woW(1:nefpvals,ii) = vel.V2woW(gd);

	% Compute and save vertical velocities and piston pos
	% Once finished, water velocity can be computed from either wp-wr or wp-wf
	
	% Constants
	g = 9.8;
	finpitch = 1; % fin pitch (pretty close to truth)
	frad = 0.126; % float radius in meters (from
	              %  diagram in 6/15/05 notes)
	%finQ = 0.605; % fin efficiency--tunable parameter (before tuning)
	if strcmp(cruise_id,'climode')|...
	      strcmp(cruise_id,'philip')|strcmp(cruise_id,'canyon'),
	  % deployments with recovery ears
	  finQ_down = 0.605*0.61; % CLIMODE tuned (with recovery loops);
	  finQ_up = 0.605*0.69; % CLIMODE tuned (with recovery loops);
	  finQ_down = 0.4770; % Philip tuned (with recovery loops);
	  finQ_up = 0.4631; % CLIMODE tuned (with recovery loops);
	elseif strcmp(cruise_id,'cblast')|strcmp(cruise_id,'eddies'),
	  % old float geometry (pvc sleeve)
	  finQ_down = 0.605*0.85; % CBLAST tuned down-profile fin efficiency
	  finQ_up = 0.605*0.98; % CBLAST tuned up-profile fin efficiency
	else,
	  % new geometry (vane ring)
	  finQ_down = 0.628;
	  finQ_up = 0.675;
	end
	if dnup(ii)==1
	  finQ = finQ_down;
	else
	  finQ = finQ_up;
	end
	
	%Cd = 0.266; % drag coefficient--tunable parameter
	Cd = 0.45; % just a guess, used ~0.5 successfully in EDDIES
	Lf = 1.233; % float length in meters (not used for drag)
	Af = 0.064; % cross-sectional area in m^2
	%Mf = float_mass(fid,'eddies');
	%Mf = float_mass(fid,cruise_id);
	[p1,T1,S1,pc1,Mf,alph1,pccpc1,Ma1] = float_ballast(fltid);

	% Now calculate velocities
	wp = vel.Wef(gd); % float pressure rate
	rotper = vel.RotP(gd); % rotation period (s)
	rhowat_ctd = sw_dens(vel.S,vel.T,vel.Pctd);

	% Drop rotper=0 and adjacent points
	bd = find(rotper<7 | [1 rotper(1:end-1)]==0 | [rotper(2:end) 1]==0);
	rotper(bd) = NaN*bd;
	%min(rotper)
	
        if dnup(ii)==1				     
	  [Pctd,j] = sort(vel.Pctd);
	else
	  [Pctd,j] = sort(-vel.Pctd);
	  Pctd = -Pctd;
	end
	[Pctd,i] = unique(Pctd);
	% So j are the sorted indices of vel.Pctd (i.e., increasing for a
	% downcast and decreasing for an upcast); i are the indices of the j
	% matrix which remove duplicate P values. Of course if the CTD
	% measurements at duplicate P values are really different, this
	% method will select the first one rather than an average, which
	% might be better.
	ji = j(i);
	rhowat_ctd = rhowat_ctd(ji);
	Tctd = vel.T(ji);
	
        % interpolate to EF grid		   
	warning off;
	rhowat = interp1(Pctd,rhowat_ctd,vel.Pef(gd),'linear');
	Twat = interp1(Pctd,Tctd,vel.Pef(gd),'linear');
	warning on;
	
	% Model
	%pc = efp.PISTON_C0; % piston position
	%pc = [NaN efp.PISTON_C0(1:end-1)]; % piston position
                            % coefficient 0 (mean)
	if isfield(vel,'pc_efp'),
	  pc = [NaN vel.pc_efp(gd(1:end-1))];
	
	  %rhof = float_den(fid,Twat,vel.Pef,pc,'eddies');
	  rhof = float_den(fltid,Twat,vel.Pef(gd),pc);
	  drho = rhowat - rhof;
	  wf = sign(drho).*sqrt((g*Mf*(1e-3)*abs(drho))./(Cd*Af*rhowat.*rhof));
        else
	  [wf,pc] = deal(vel.Pef(gd).*NaN);
	end

	% Note that John's processed files don't have any indicator of
	% whether the float is rotating clockwise or counterclockwise. For
	% now, assume rotation is always clockwise on down profiles and
	% counterclockwise on up. Eliminates the identification of reversals
	% due to strong vertical velocity, but these are pretty unlikely
	% anyway...
	
	if mod(hprof,2)==1
	  rdir = -1;
	else
	  rdir = 1;
	end
	wr = rdir*(2*pi*frad)./(rotper*finpitch*finQ); % velocity inferred
                                                  % from rotation rate

	% and add a few more variables to keep track of pressure rate and
	% piston position on the CTD grid (which seems to save more often
	% during the hold phase)
	%
	% Question: why are the piston values on the CTD and EFP grids offset
	% in time by 20-40s (and varying with depth)? and which one is
	% correct? Could the EFP be reporting the final piston position (at the end
	% of the 50s recording interval) along with the midpoint time? This
	% makes pc timeseries match much better, but EFP still has
	% variations. Could be the correct averages on a shorter time
	% grid. In fact, most accurate timeseries might be to combine the
	% two (or even use higher-order poly parameters from efp?). Then CTD
	% diff timegrid could be a bin average or interpolation from this timeseries.
        dpdt_ctd = -diff(vel.Pctd)./diff(vel.ctd_mlt)/86400;
        ctd_tav = (vel.ctd_mlt(1:end-1)+vel.ctd_mlt(2:end))/2;
	%pca_av = [(ctd.pca(2:end-1)+ctd.pca(3:end))/2 ctd.pca(end)];
	%pca_av = [(vel.pca(2:end-1)+vel.pca(3:end))/2 vel.pca(end)];
	%pca_av = [(vel.pca(1:end-1)+vel.pca(2:end))/2];
	if isfield(vel,'pc_ctd'),
	  pca_av = [(vel.pc_ctd(1:end-1)+vel.pc_ctd(2:end))/2];
	else
	  pca_av = vel.Pctd(1:end-1)*NaN;
	end
	p_av = (vel.Pctd(1:end-1)+vel.Pctd(2:end))/2;
	
	p_av(lints) = NaN*lints;
	pca_av(lints) = NaN*lints;
	ctd_tav(lints) = NaN*lints;
	dpdt_ctd(lints) = NaN*lints;
	
	if do_wplots
	  figure(1)
	  % plots vs. pressure
	  subplot(311)
	  %plot(vel.Pef,pc,'r.-',vel.Pctd,ctd.pca,'b.-');
	  if isfield(vel,'pc_ctd'),
	    plot(vel.Pef(gd),pc,'r.-',vel.Pctd,vel.pc_ctd,'b.-');
	  else
	    plot(vel.Pef(gd),pc,'r.-');
	  end
	  grid on
	  ylabel('piston counts');
	  title(['EM-APEX ' num2str(fid) '; Profile #' num2str(hprof)]);
		       
	  subplot(312)
	  plot(vel.Pef(gd),wr,'m.-',p_av,dpdt_ctd,'b.-',vel.Pef(gd),wp,'r',vel.Pef(gd),wf,'g');
	  grid on
	  ylabel('float w (m/s)');
	  legend('w_{rotper}','w_{dpdt-CTD}','w_{dpdt-EF}','w_{float\Delta\rho}',4);

	  subplot(313)
	  plot(vel.Pef(gd),wp-wf,vel.Pef(gd),wp-wr,'r');
	  grid on
	  xlabel('pressure (dbar)')
	  ylabel('water w (m/s)');
	  legend('float model','current meter',4);
		   
	  figure(2)
	  % plots vs. time
	  % Note: could also include remaining curves in fig.1 including
	  % w_float_Delta_rho and "water w" panel. Then this fig would be a
	  % 4x1 subplot.
	  
	  subplot(311)
	  plot(vel.ctd_mlt,vel.Pctd,'c-',ctd_tav,p_av,'b.-',vel.efp_mlt,vel.Pef,'r.-');
	  grid on; axis ij;
	  datetick; ylabel('pressure');
	  title('Red=EF grid; Blue=CTD grid');
	  
	  subplot(312)
	  %plot(vel.efp_mlt,pc,'r.-',vel.ctd_mlt,ctd.pca,'c-',ctd_tav,pca_av,'b.-'); grid on;
          if isfield(vel,'pc_ctd'),
	  plot(vel.efp_mlt(gd),pc,'r.-',vel.ctd_mlt,vel.pc_ctd,'c.-',ctd_tav,pca_av,'b.-');
	  else
	  plot(vel.efp_mlt(gd),pc,'r.-',ctd_tav,pca_av,'b.-');
	  end
	  grid on;
	  datetick; ylabel('piston counts');
	  legend('EF grid','CTD grid','pca_av');
	  
	  subplot(313)
	  %plot(vel.efp_mlt,wp,'r.-',ctd_tav,dpdt_ctd,'b.-');
	  plot(vel.efp_mlt(gd),wp,'r.-',ctd_tav,dpdt_ctd,'b.-',vel.efp_mlt(gd),wr,'m.-');
	  grid on;
	  datetick; ylabel('dpdt');
	  legend('w_{dpdt-EF}','w_{dpdt-CTD}','w_{rotper}');
	  
	  disp('paused');pause
	end % 	if do_wplots
	
	Wp(1:nefpvals,ii) = wp;
	Wr(1:nefpvals,ii) = wr;
	Wf(1:nefpvals,ii) = wf;
	ppos(1:nefpvals,ii) = pc;
	
	% Compute variables on CTD midpoint (dpdt) grid
	Wp_ca(1:npvals-1,ii) = dpdt_ctd;
	ppos_ca(1:npvals-1,ii) = pca_av;
	P_ca(1:npvals-1,ii) = p_av;
	mlt_ca(1:npvals-1,ii) = ctd_tav;
      end % if ~ctdonly


      % Use shallowest 3 UTC and P values to extrapolate to the surfacing
      % time
      npoints = 3; % number of pressure values to fit to
      %pthresh = 15; % only include profiles that reach 15m depth
      pord = 1; % polynomial fit order

      if dnup(ii)==2,
	if npvals<3
	  error(['Dead end (821): Only ' num2str(npvals) ' pressure values in profile ' ...
		num2str(hprof) '.']);
	end % if npvals<3
	pfit = vel.Pctd(end-npoints+1:end);
	timfit = vel.ctd_mlt(end-npoints+1:end);
	c = polyfit(pfit,timfit,pord);
      
%	  if min(pfit)<pthresh,
	    % only extrapolate if profile got within 15m of the surface... this
	    % mainly distinguishes between real up profiles and yoyos
	  if surfaced(ii),
	    % if profile made it to the surface
	    %utc_up(ii) = c(2);
	    utc_up(ii) = polyval(c,0);
	  end

	if isnan(utc_desc)&hprof>2,
	  % Here compute synthetic down profile for Park (Argo) Mode profiles (up
	  % without previous down), since utc_desc should only be NaN here directly
	  % after an up that reached the surface.
	  if verbose
	    disp(['Starting subsurface interval at up profile ' num2str(hpid(ii))...
		  ' using estimated descent after previous up.' ]);
	  end
	  %utc_desc = utc_gps(ii-1)+20*60/86400; % guess: 20 min after
          %                                      % first good GPS fix.
	  utc_desc = utc_up(ii-1)+21*60/86400; % guess: 21 min after
                                                % previous surfacing.
	  % (Estimated on the basis of CBLAST performance.)
	  % Note that the definition based on utc_up is better if GPS is
	  % missing, since it depends only on the APF9 clock.
	  %
	  % To refine this, need to know what tasks were included in the
	  % surface interval (i.e., second GPS?) and how much data was
	  % sent (or really, just the Telemetry timeout). Do phone logs
	  % contain info on connection time?
	  %
	  % Record this as down time as well, since there's no actual down
	  % profile. This suggests a way to identify "Argo mode" profiling, since
	  % these are the only profiles that include both utc_down and utc_up values.
	  utc_down(ii) = utc_desc;			    
					    
      % Also need to compute pressure timeseries and interpolated up velocity
      % profile from estimated park profiling behavior
      
    end % if isnan(utc_desc)&ii>1
    
    elseif dnup(ii)==1,
      % compute utc_down for all down profiles, even if never at surface
      % No, scratch that. Better to leave it NaN unless previous profile
      % actually made it to the surface.
      if npvals<3
	error(['Dead end (872): Only ' num2str(npvals) ' pressure values in profile ' ...
	      num2str(hprof) '.']);
      end
      pfit = vel.Pctd(1:npoints);
      timfit = vel.ctd_mlt(1:npoints);
      c = polyfit(pfit,timfit,pord);
      %utc_desc(ii) = c(2);
      if surfaced(ii),
	% Profile came from surface
	utc_down(ii) = polyval(c,0);
      end
      if isnan(utc_desc)
	% Save the last real departure from the surface on first down after
	% surface interval (when utc_desc was reset)
	%
	% Now need to decide when to compute descent position. Also, this
	% is going to be tougher with the later "Argo" profiles which have
	% no down data. Need to make some guesses... For up profiles with
	% no down, probably should save
	% an estimated down time and position anyway. Also need to make up
	% a velocity timeseries based on the up profile and the
	% predicted Park mode pressure timeseries based on the piston position and
	% other parameters
	% in effect.
	if verbose
	  disp(['Starting subsurface interval at profile ' num2str(hpid(ii))]);
	end
	if hprof==1,
	  utc_desc = ddep;
	  disp(['Float ' fltid ' took ' num2str((utc_down(ii)-utc_desc)*24*60) ...
		' min longer to sink than expected from CTD pressure fit.']);
	else
	  utc_desc = utc_down(ii);
	end
      end

    else
	error('dnup must be 1 or 2!');
      end % if dnup(ii)==2,

      if isnan(utc_desc_vbs),
	if verbose
	  disp(['Starting vbarstar interval at profile ' num2str(hpid(ii))]);
	end
	utc_desc_vbs = utc_desc;
      end
      
      if isnan(lat_desc_vbs)&isfinite(utc_desc_vbs)&exist('cx'),
	% extrapolate surface drift forward to utc_desc_vbs
	x_desc = polyval(cx,(utc_desc-t0)*86400);
	y_desc = polyval(cy,(utc_desc-t0)*86400);
	lon_desc = x_desc/lat2m/cos(lat0*pi/180) + lon0;
	lat_desc = y_desc/lat2m + lat0;

	lat_down(ii) = lat_desc;
	lon_down(ii) = lon_desc;
	lat_desc_vbs = lat_desc;
	lon_desc_vbs = lon_desc;
	
	if do_vbsplots
	  figure(3); clf
	  plot(ogps.lon(oi1),ogps.lat(oi1),'c.-',...
	      lon_desc,lat_desc,'ms');
	  grid on; hold on
	  set(gca,'DataAspectRatio',[cos(lat0*pi/180) 1 1]);
	  xlabel('lon');ylabel('lat');
	  title(['ms: descent position; c.:previous gps']);
	  %disp('paused with first part of vbs plot 3'); pause
	end
	
	clear cx cy t0 lat0 lon0 ogps oi1
      elseif hprof==1 & ~isnan(lat_deploy(fltidind)),
	% for first profile, use deployment position
	lat_desc = lat_deploy(fltidind);
	lon_desc = lon_deploy(fltidind);
	lat_desc_vbs = lat_desc;
	lon_desc_vbs = lon_desc;
      end
      
      nvmin = 2;
%      if ~ctdonly & nvel(ii)>2,
      if ~ctdonly & nvel(ii)>nvmin,
	% Accumulate velocity timeseries for absolute velocity (needs to be
	% after utc_desc is calculated for time base)
	% if dnup(ii)==2 & (ii==1 | dnup(ii-1)==2),
	
	clear Tg Pg Ugap Vgap Tshal Pshal Ushal Vshal
	
	if argo_mode(ii)
	  % For up profiles without a previous down ('Argo mode'), construct
	  % a synthetic down profile from expected pressure record and up
	  % velocity data.

	  if 0
	  % First guess: fake it and just assume a linear descent from
	  % surface at utc_desc to ema.Pef(1) at ema.efp_mlt(1)
	  Pdesc = [0:2:vel.Pef(1)];
	  Tdesc = interp1([0 vel.Pef(1)],[utc_desc vel.efp_mlt(1)],...
	      Pdesc,'linear');
	  Timacc = [Timacc; Tdesc(2:end-1)'-utc_desc];
	  end
	
	  if 0
	  % Next guess: Assume linear descent at 0.08 m/s to max depth (500m),
	  % then park
	  dspd = 0.08;
	  Pdesc = [0:2:vel.Pef(1)];
	  Tdesc = Pdesc/dspd/86400;
	  Timacc = [Timacc; Tdesc(2:end-1)'];
	  end

	  % Third guess: Assume exponential descent with 0.2 day decay scale
	  % to 500m
	  tau = 0.2;
	  gd = find(isfinite(vel.Pef)&isfinite(vel.efp_mlt));
	  Tdesc = [0:100/86400:vel.efp_mlt(gd(1))-utc_desc-100/86400];
	  Pdesc = vel.Pef(gd(1))*(1 - exp(-Tdesc/tau));
	  
	  gd2 = find(isfinite(Ug)&isfinite(Vg));
	  Udesc = interp1(vel.Pef(gd(gd2)),Ug(gd2),Pdesc,'linear');
	  Vdesc = interp1(vel.Pef(gd(gd2)),Vg(gd2),Pdesc,'linear');

	  Timacc = [Timacc; Tdesc(2:end-1)'];
	  Pacc = [Pacc; Pdesc(2:end-1)'];
	  Uacc = [Uacc; Udesc(2:end-1)'];
	  Vacc = [Vacc; Vdesc(2:end-1)'];
	  %Timacc = [Tdesc(2:end-1)'];
	  %Pacc = [Pdesc(2:end-1)'];
	  %Uacc = [Udesc(2:end-1)'];
	  %Vacc = [Vdesc(2:end-1)'];
	end % 	if argo_mode(ii)

	% Set up initial U,V,P,t to accumulate (subject to modification by
	% patching in upper ocean or hold-depth gap)
	gd = find(isfinite(Pef(:,ii)) & isfinite(UTCef(:,ii)));
	Uga = U(gd,ii);
	Vga = V(gd,ii);
	Ta = UTCef(gd,ii) - utc_desc_vbs;
        Pa = Pef(gd,ii);

	bd = find(isnan(Ta)|isnan(Pa));
	Pa(bd) = []; Vga(bd) = []; Uga(bd)=[]; Ta(bd)=[];

	if dnup(ii)==1 & gotuni,
	  % Near-surface substitution. First of all, only enter this on down
	  % profiles and only the second time around (gotuni==1)

	  % If down profile starts much deeper than up, extrapolate from 
	  % shallowest measurement using up-profile shear
	  
	iup = ii+1; % default to next up profile (although previous up is usually
	% closer in time)

	if ii>1 & flid(ii-1)==flid(ii) & depl(ii-1)==depl(ii),
	  % for all but first profile of mission, use previous up
	  iup = ii-1;
	end
	% but check to make sure velocity has at least 2 good values; maybe
	% even iterate here?
	igd = find(isfinite(uni.U(:,iup)) & uni.U(:,iup)~=0);
	if length(igd)<3 & iup==ii-1 & length(uni.hpid)>ii;
	  iup=ii+1;
	  igd = find(isfinite(uni.U(:,iup)) & uni.U(:,iup)~=0);
	end
	
	if length(uni.hpid)>=max(iup,ii) & ...
	      min(uni.Pef(:,ii))-min(uni.Pef(:,iup))>10 & dnup(iup)==2,
		      
	  [mPd,imPd] = min(uni.Pef(:,ii)); % min P on down
	  ishal = find(uni.Pef(:,iup) < mPd & uni.Pef(:,iup)~=0);
	  u0 = interp1(uni.Pef(igd,iup),uni.U(igd,iup),mPd,'linear');
	  v0 = interp1(uni.Pef(igd,iup),uni.V(igd,iup),mPd,'linear');

	  
	  Ushal = flipud(uni.U(ishal,iup)-u0+uni.U(imPd,ii));
	  Uga = [Ushal; Uga];
	  Vshal = flipud(uni.V(ishal,iup)-v0+uni.V(imPd,ii));
	  Vga = [Vshal; Vga];
	  %Pa = [flipud(uni.Pef(ishal,iup)); vel.Pef(gd)'];
	  Pshal = flipud(uni.Pef(ishal,iup));
	  Pa = [Pshal; Pa];

%	  Tshal = interp1([0 mPd],...
%	      [utc_desc min(vel.efp_mlt)]-utc_desc_vbs,...
%	      flipud(uni.Pef(ishal,iup)),'linear');
          if surfaced(ii)
	    % assume that if surfaced(ii) then the vbs interval is starting
	    % with the current down profile so get time by interpolating
	    % inserted pressures from the surface.
	    Tshal = interp1([0 mPd],...
		[utc_desc uni.UTCef(imPd,ii)]-utc_desc_vbs,...
		flipud(uni.Pef(ishal,iup)),'linear');
	  else
	    % if not coming from surface, interpolate time corresponding to
	    % inserted pressures from the
	    % shallowest (last) t,P point of the previous (up) profile
	    [plast,iplast] = min(uni.Pef(:,ii-1));
	    tlast = uni.UTCef(iplast,ii-1);
	    Tshal = interp1([plast mPd],...
		[tlast uni.UTCef(imPd,ii)]-utc_desc_vbs,...
		flipud(uni.Pef(ishal,iup)),'linear');
	  end
	  %Ta = [Tshal; vel.efp_mlt(gd)' - utc_desc_vbs];
	  Ta = [Tshal; Ta];

	  
	  %Uacc = [Uacc; Uga];
	  %Vacc = [Vacc; Vga];
	  %Timacc = [Timacc; Ta];
	  %Pacc = [Pacc; Pa];

	  elseif length(uni.hpid)>=max(iup,ii) & ...
		length(igd)>3 & isempty(find(isfinite(uni.Pef(:,ii)))) & dnup(iup)==2,
	    % For profiles with CTD data but no EF at all, insert EF from iup
	    % (ideally the nearest good up prof) profile (tuned for Philex
	    % 3305c hp351... works for any others?). Since CTD data includes
	    % best info about float trajectory, keep time and depth grids and
	    % interpolate EF to those. (Maybe even better would be to
	    % interpolate to a finer time grid first?

	    disp([fltid ': Inserting hp' num2str(hpid(iup)) ' velocity into missing hp' num2str(hpid(ii))]);
	    
	    gd = find(isfinite(UTC(:,ii)));
	    Ta = UTC(gd,ii)-utc_desc_vbs;
	    Pa = P(gd,ii);
	    Uga = interp1(uni.Pef(igd,iup),uni.U(igd,iup),Pa,'linear');
	    Vga = interp1(uni.Pef(igd,iup),uni.V(igd,iup),Pa,'linear');
	    
	    [mPd,imPd] = min(P(gd,ii)); % min P on down
	  ishal = find(uni.Pef(:,iup) < mPd & uni.Pef(:,iup)~=0);

	  Ushal = flipud(uni.U(ishal,iup));
	  Uga = [Ushal; Uga];
	  Vshal = flipud(uni.V(ishal,iup));
	  Vga = [Vshal; Vga];
	  %Pa = [flipud(uni.Pef(ishal,iup)); vel.Pef(gd)'];
	  Pshal = flipud(uni.Pef(ishal,iup));
	  Pa = [Pshal; Pa];

	    Tshal = interp1([0 mPd],...
		[utc_desc UTC(gd(imPd),ii)]-utc_desc_vbs,...
		Pshal,'linear');
	    Ta = [Tshal; Ta];
	    
	end % 	if length(uni.hpid)>=ii+1 & min(uni.Pef(:,ii))-min(uni.Pef(:,ii+1))>10,
      end % if if dnup(ii)==1 & gotuni
	
	bd = find(isnan(Ta)|isnan(Pa));
	Pa(bd) = []; Vga(bd) = []; Uga(bd)=[]; Ta(bd)=[];

	if dnup(ii)==1 & gotuni & max(diff(UTC(:,ii))*86400)>290,
	  % If this is a down-profile with a "hold" interval, need to add in
	  % CTD pressure and time information along with a synthesized
	  % velocity profile (probably using previous up prof). This is
	  % critical, since the great majority of the float's lifetime may be
	  % spent in these holds so this may be the dominant contribution to
	  % the GPS velocity!
	  %pgd = find(isfinite(vel.Pef) & isfinite(vel.efp_mlt));
%	  pgd = find(isfinite(Pef(:,ii)) & isfinite(UTCef(:,ii)));
	  cgd = find(isfinite(P(:,ii)) & isfinite(UTC(:,ii)));
	  
	  % find the biggest gap in vel.efp_mlt and fill it with CTD time and
	  % pressure data, then interpolate pressures to previous up profile
	  % velocity (or maybe an average of previous and following
	  % up-velocity?)
	  
	  Ta2 = [UTC(cgd(1),ii)-utc_desc_vbs; Ta; UTC(cgd(end),ii)-utc_desc_vbs];
	  %[dp,ig] = max(diff(vel.efp_mlt(pgd)));
%	  [dp,ig] = max(diff(UTCef(pgd,ii)));
	  [dp,ig] = max(diff(Ta2));
	  
%	  icg = find(vel.ctd_mlt>vel.efp_mlt(pgd(ig)) & ...
%	      vel.ctd_mlt<vel.efp_mlt(pgd(ig+1)));
	  icg = find(UTC(cgd,ii)>Ta2(ig)+utc_desc_vbs & ...
	      UTC(cgd,ii)<Ta2(ig+1)+utc_desc_vbs);
	  
	  if ~isempty(icg),
	  Tg = UTC(cgd(icg),ii);
	  Pg = P(cgd(icg),ii);
	  iiprev = max(find(uni.flid==fid & uni.depl==depl(ii) & ...
	      mod(uni.hpid,2)==0 & uni.nvel>nvmin & uni.hpid<hpid(ii) & max(uni.P)>=max(Pg)));
	  iinext = min(find(uni.flid==fid & uni.depl==depl(ii) & ...
	      mod(uni.hpid,2)==0 & uni.nvel>nvmin & uni.hpid>hpid(ii) & max(uni.P)>=max(Pg)));
	  if isempty(iiprev)&isempty(iinext),
	    % This is the deepest profile of the deployment AND there is no
	    % bottom segment... maybe only true for Philex:3304c?
	    % Use the deepest up profile available
	    mP = max(uni.P);
	    [m2P,im2P] = max(mP(find(uni.flid==fid & uni.depl==depl(ii) & ...
		mod(uni.hpid,2)==0 & uni.hpid~=hpid(ii))));
	    iiprev = max(find(uni.flid==fid & uni.depl==depl(ii) & ...
		mod(uni.hpid,2)==0 & uni.hpid<hpid(ii) & max(uni.P)>=m2P));
	    iinext = min(find(uni.flid==fid & uni.depl==depl(ii) & ...
		mod(uni.hpid,2)==0 & uni.hpid>hpid(ii) & max(uni.P)>=m2P));
	    gd3 = find(Pg<=m2P);
	    Pg = Pg(gd3); Tg = Tg(gd3);
	  end
	  
	  if isempty(iiprev),
	    gd2 = find(isfinite(uni.Pef(:,iinext))&isfinite(uni.U(:,iinext)));
	    Ugap = interp1(uni.Pef(gd2,iinext),uni.U(gd2,iinext),Pg,'linear');
	    Vgap = interp1(uni.Pef(gd2,iinext),uni.V(gd2,iinext),Pg,'linear');
	  elseif isempty(iinext),
	    gd1 = find(isfinite(uni.Pef(:,iiprev))&isfinite(uni.U(:,iiprev)));
	    Ugap = interp1(uni.Pef(gd1,iiprev),uni.U(gd1,iiprev),Pg,'linear');
	    Vgap = interp1(uni.Pef(gd1,iiprev),uni.V(gd1,iiprev),Pg,'linear');
	  else
	    gd1 = find(isfinite(uni.Pef(:,iiprev))&isfinite(uni.U(:,iiprev)));
	    gd2 = find(isfinite(uni.Pef(:,iinext))&isfinite(uni.U(:,iinext)));
	    
	    if length(gd1)<2
	      Ug2 = interp1(uni.Pef(gd2,iinext),uni.U(gd2,iinext),Pg,'linear');
	      Vg2 = interp1(uni.Pef(gd2,iinext),uni.V(gd2,iinext),Pg,'linear');
	      Ug1 = Ug2; Vg1 = Vg2; w1 = 0; w2 = 1;
	    elseif length(gd2)<2
	      Ug1 = interp1(uni.Pef(gd1,iiprev),uni.U(gd1,iiprev),Pg,'linear');
	      Vg1 = interp1(uni.Pef(gd1,iiprev),uni.V(gd1,iiprev),Pg,'linear');
	      Ug2 = Ug1; Vg2 = Vg1; w2 = 0; w1 = 1;
	    else
	      Ug1 = interp1(uni.Pef(gd1,iiprev),uni.U(gd1,iiprev),Pg,'linear');
	      Vg1 = interp1(uni.Pef(gd1,iiprev),uni.V(gd1,iiprev),Pg,'linear');
	      Ug2 = interp1(uni.Pef(gd2,iinext),uni.U(gd2,iinext),Pg,'linear');
	      Vg2 = interp1(uni.Pef(gd2,iinext),uni.V(gd2,iinext),Pg,'linear');

	      t1 = mymean(uni.UTCef(:,iiprev));
	      t2 = mymean(uni.UTCef(:,iinext));
	      w1 = (t2-Tg)/(t2-t1);
	      w2 = (Tg-t1)/(t2-t1);
	    end
	    
	    Ugap = Ug1.*w1+Ug2.*w2;
	    Vgap = Vg1.*w1+Vg2.*w2;
	  end	    
	  
	  %Ta = [vel.efp_mlt(pgd(1:ig))'; Tg'; vel.efp_mlt(pgd(ig+1:end))'] - utc_desc_vbs;
%	  Ta = [UTCef(pgd(1:ig),ii); Tg; UTCef(pgd(ig+1:end),ii)] - utc_desc_vbs;
	  %Pa = [vel.Pef(pgd(1:ig))'; Pg'; vel.Pef(pgd(ig+1:end))'];
%	  Pa = [Pef(pgd(1:ig),ii); Pg; Pef(pgd(ig+1:end),ii)];
%	  Uga = [Ug(pgd(1:ig)); Ugap; Ug(pgd(ig+1:end))];
%	  Vga = [Vg(pgd(1:ig)); Vgap; Vg(pgd(ig+1:end))];

          if ig==1,
	    Ta = [Tg - utc_desc_vbs; Ta(ig:end)];
	    Pa = [Pg; Pa(ig:end)];
	    Uga = [Ugap; Uga(ig:end)];
	    Vga = [Vgap; Vga(ig:end)];
          elseif ig==length(Ta)+1,
	    Ta = [Ta(1:ig-1); Tg - utc_desc_vbs];
	    Pa = [Pa(1:ig-1); Pg];
	    Uga = [Uga(1:ig-1); Ugap];
	    Vga = [Vga(1:ig-1); Vgap];
	  else
	    Ta = [Ta(1:ig-1); Tg - utc_desc_vbs; Ta(ig:end)];
	    Pa = [Pa(1:ig-1); Pg; Pa(ig:end)];
	    Uga = [Uga(1:ig-1); Ugap; Uga(ig:end)];
	    Vga = [Vga(1:ig-1); Vgap; Vga(ig:end)];
	  end

%	  Uacc = [Uacc; Uga];
%	  Vacc = [Vacc; Vga];
%	  Timacc = [Timacc; Ta];
%	  Pacc = [Pacc; Pa];

	  % So far, looks like this procedure works well for long holds, but
	  % not so good for short holds. Why? Not clear, but somehow the
	  % profile-to-profile variability makes interpolation between before
	  % and after worse than just sticking with original profile
	  % (especially when down segments overlap instead of having a big
	  % gap). Maybe we're just talking about an inevitable 2 cm/s noise level
	  % or so because of high-frequency variability?
	  %
	  % In the end, is the resulting ubs,vbs estimate more like an
	  % average over the entire interval or just an average of samples at
	  % the beginning and end? Not quite either, since GPS vel is a full
	  % avg while EF vel is 2 samples! Resulting VBS is a hybrid and
	  % accuracy depends on high-freq noise as well as long-time
	  % averaging. So the next question is how do we evaluate ubs,vbs
	  % noise level? From observed variability and short-long term differences?
	  
	  end % if ~isempty(icg)
	end % 	if dnup(ii)==1 & gotuni & max(diff(UTC(:,ii))*86400)>290,

	if do_plots_each
	    figure(4); clf
	    subplot(221)
	    plot(Timacc+utc_desc_vbs,Pacc,'b.',Ta+utc_desc_vbs,Pa,'r.-');
	    if exist('Tg'),
	      hold on; plot(Tg,Pg,'m.-');
	    end
	    if exist('Tshal'),
	      hold on; plot(Tshal+utc_desc_vbs,Pshal,'mx-');
	    end
	    grid on; axis ij; datetick
	    subplot(223)
	    plot(Timacc+utc_desc_vbs,Uacc,'r.-',Ta+utc_desc_vbs,Uga,'r+-');
	    if exist('Tg'),
	      hold on; plot(Tg,Ugap,'m.-');
	    end
	    if exist('Tshal'),
	      hold on; plot(Tshal+utc_desc_vbs,Ushal,'mx-');
	    end
	    hold on
	    plot(Timacc+utc_desc_vbs,Vacc,'b.-',Ta+utc_desc_vbs,Vga,'b+-');
	    if exist('Tg'),
	      hold on; plot(Tg,Vgap,'c.-');
	    end
	    if exist('Tshal'),
	      hold on; plot(Tshal+utc_desc_vbs,Vshal,'cx-');
	    end
	    grid on; datetick
	    subplot(122)
	    plot(Uacc,Pacc,'r.-',Uga,Pa,'r+-');
	    hold on;
	    plot(Vacc,Pacc,'b.-',Vga,Pa,'b+-');
	    title(['EM-APEX ' num2str(fltid) '; hpid:' num2str(hprof)...
	    '; ssid:' num2str(ssid(ii))]);
	    if exist('Tg'),
	    if ~isempty(iiprev)
	      plot(uni.U(gd1,iiprev),uni.Pef(gd1,iiprev),'r-.',uni.V(gd1,iiprev),uni.Pef(gd1,iiprev),'b-.');
	    end
	    if ~isempty(iinext)
	      plot(uni.U(gd2,iinext),uni.Pef(gd2,iinext),'r:',uni.V(gd2,iinext),uni.Pef(gd2,iinext),'b:');
	    end
	      hold on; plot(Ugap,Pg,'m.-',Vgap,Pg,'c.-');
	    end
	    if exist('Tshal'),
	      hold on; plot(Ushal,Pshal,'mx-',Vshal,Pshal,'cx-');
	    end
	      axis ij; grid on;
	    %pause
	    
	    figure(5)
	    % one more plot, this time showing E1 and E2 comparison
	    clf
	    plot(U1(:,ii),Pef(:,ii),'r.-',U2(:,ii),Pef(:,ii),'m.-');
	    hold on;
	    plot(V1(:,ii),Pef(:,ii),'b.-',V2(:,ii),Pef(:,ii),'c.-');
	    title(['EM-APEX ' num2str(fltid) '; hpid:' num2str(hprof)...
		  '; ssid:' num2str(ssid(ii)) ' (r,m: U1,U2; b,c: V1,V2)']);
	    axis ij; grid on;

	  end

	bd = find(isnan(Ta)|isnan(Pa));
	Pa(bd) = []; Vga(bd) = []; Uga(bd)=[]; Ta(bd)=[];

	Uacc = [Uacc; Uga];
	Vacc = [Vacc; Vga];
	Timacc = [Timacc; Ta];
	Pacc = [Pacc; Pa];

      end % if ~ctdonly & nvel(ii)>2,


      if length(Uacc)~=length(Timacc),
	error('Uacc and Timacc don''t match!');
      end
     
      utc_beg(ii) = vel.ctd_mlt(1);
      
      utc_dep(ii) = ddep;
      utc_rec(ii) = drec;
      
      % Save VIT data
      
      % Save GPS data (indexed by hpid, so ...
      % See lines 404-477 in unifiles_fran.m

%      if strcmp(cruise_id,'cblast'),
%	gotgps = vel.gotgps;
%      elseif ~isempty(find(gps.hpid==hprof)),
      if ~isempty(find(gps.hpid==hprof)),
	% Not sure if this is exactly how the cblast files determined the
	% gotgps variable, but probably close.
	gotgps = 1;
      else
	gotgps = 0;
      end

      if verbose
      disp([num2str(ii) ' ' num2str(fltid) ' hpid:' num2str(hprof)...
	    ' ssid:' num2str(profnum) ' gotgps:' num2str(gotgps)]);
      end

      i1 = [];
      if gotgps
	% GPS file exists (i.e., we made it to the surface) but is there any good data?

%	if strcmp(cruise_id,'cblast'),
	if 0
	  % Use "raw" GPS file instead of GPS field in "processed" file. Why?
	  % (1) GPS files contain APF9/GPS time comparison and (2) GPS files
	  % sometimes contain a few additional fixes. Only drawback is that
	  % GPS files have incorrect hpid once Park mode profiling
	  % starts. Use lastfastprofile to correct numbers after this.
	  %if fid==1634, % since only 1634 contains APF9/GPS comparison

	  lfp = lastfastprof(fltidind);
	  if hprof>lfp,
	    rprof = (hprof-lfp)/2+lfp;
	  else
	    rprof = hprof;
	  end
	  gpsname = ['ema-' rfltid(1:4) '-gps-' sprintf('%04d',rprof) '.mat'];
	  %gpsname = ['ema-' fltid(1:4) '-gps-' sprintf('%04d',hprof) '.mat']
	  GPS = load([rawdir gpsname]);
	  %apf9_terr(ii) = median((gps.mlt_apf9-gps.mlt_gga)*86400);
	
	  if do_navplots
	    clf
	    plot(GPS.lon,GPS.lat,'m+',ema.GPS.lon,ema.GPS.lat,'o'); grid on;
	    title('+: GPS raw file; o: ema proc file');
	    disp([num2str(fltid) ' hpid:' num2str(hprof) ' ema start: '...
		  datestr(GPS.utc(1)) ' gps start: ' datestr(GPS.mlt_gga(1))]);
	    disp('showing all GPS data; paused'); pause
	  end
	  
	  %end
	else

	  if do_navplots
	    % plot full GPS track for this float
	    figure(1); clf
	    gd = find(gps.stat==1&gps.nsat>3);
	    plot(gps.lon(gd),gps.lat(gd),'m+'); grid on;
	    title('+: GPS file');
	    disp([num2str(fltid) ' hpid:' num2str(hprof) ' APF9 start: '...
		  datestr(gps.mlt_apf9(1)) ' GGA start: ' datestr(gps.mlt_gga(1))]);
	    disp('showing all GPS data; paused'); pause
	  end
	  
	end % if strcmp(cruise_id,'cblast'),

	profnum = profnum + 1;
	i1 = find(gps.hpid==hprof&gps.stat==1&gps.nsat>3);
	apf9_terr(ii) = median((gps.mlt_apf9(i1)-gps.mlt_gga(i1))*86400);
      
      end % if gotgps

      
      % other possible criteria: hdop>2.5
      if ~isempty(i1),
	% GPS file exists and has good data... Save position, drift, and
	% APF9/GPS difference if available;
	% compute vbarstar from accumulated velocities and displacement from
	% last surfacing.

	% Note: Could subtract one hour from GPS time here. Was needed for CBLAST
	% and EDDIES, but as of 1/22/09 doesn't seem to be.
	% utc_gps(ii) = gps.mlt_gga(i1(1))-1/24;

	utc_gps(ii) = gps.mlt_gga(i1(1));
	lat_gps(ii) = gps.lat(i1(1));
	lon_gps(ii) = gps.lon(i1(1));

	if 0
	if ~gotuni | isnan(uni.botdep(ii)),
	disp(['interpolating bathymetry for float ' num2str(fid) '...']);
	botdep(ii) = interp2(tlon,tlat,-topo,lon_gps(ii),lat_gps(ii));
	disp(['...done']);
        else 
	  botdep(ii) = uni.botdep(ii);
        end
        end
	
	median_gpstime(ii) = median(gps.mlt_gga(i1));
	first_apf9time(ii) = gps.mlt_apf9(i1(1));
	median_apf9time(ii) = median(gps.mlt_apf9(i1));
      
      	if do_navplots
	  % highlight latest profile on full track (fig1) and add separate
	  % figure (fig2) with just current surface fixes (before and after
	  % phone call)
	  figure(1); hold on; plot(gps.lon(i1),gps.lat(i1),'b.',...
	      lon_gps(ii),lat_gps(ii),'rs','markersize',10); hold off;
	  
	  figure(2); clf
	  hold on; plot(gps.lon(i1),gps.lat(i1),'b.',...
	      lon_gps(ii),lat_gps(ii),'rs','markersize',10); hold off;
	  xlabel('lon');ylabel('lat');title('s: saved position; dots: good fixes');
	  grid on
	  disp(['clock error: ' num2str(apf9_terr(ii)) 's']);
	  disp(['showing ' num2str(length(i1))...
		' fixes used for surface drift (dots) '...
		'and position (square)']);
	  disp('paused'); pause
        end

	if length(i1)>=2,
	% estimate surface drift speed (use with utc_gps and lat_gps/lon_gps
	% to extrapolate back to utc_up)
	lat2m = 111.12e3;
	t0 = gps.mlt_apf9(i1(1));
	lon0 = gps.lon(i1(1));
	lat0 = gps.lat(i1(1));
	x = (gps.lon(i1) - lon0)*lat2m*cos(lat0*pi/180);
	y = (gps.lat(i1) - lat0)*lat2m; % m
	t = (gps.mlt_apf9(i1) - t0)*86400; % s
	cx = polyfit(t,x,1);
	u_sfc(ii) = cx(1);
	cy = polyfit(t,y,1);
	v_sfc(ii) = cy(1);
	gps_dist(ii) = sqrt((max(x)-min(x))^2+(max(y)-min(y))^2);
	
	x_up = polyval(cx,(utc_up(ii)-t0)*86400);
	y_up = polyval(cy,(utc_up(ii)-t0)*86400);
	lon_up(ii) = x_up/lat2m/cos(lat0*pi/180) + lon0;
	lat_up(ii) = y_up/lat2m + lat0;
	
	ogps = gps; oi1 = i1;
	
	% Also need to extrapolate forward to utc_desc. Keep
	% t0,lon0,lat0,cx,cy until utc_desc is known (from next descending
	% pressure record or from Argo mode punt).

	if do_vbsplots
	  figure(3); hold on;
	  plot(gps.lon(i1),gps.lat(i1),'k.-',...
	      [lon_desc lon_up(ii)],[lat_desc lat_up(ii)],'r.--');
	  plot([lon_down(dninds) lon_up(upinds)],...
	      [lat_down(dninds) lat_up(upinds)],'m.--');
	  hold off; grid on;
	  %disp('Paused with vbs plot');	  pause
	end

	if do_navplots
	  % plot distance (meters) vs. time
	subplot(211); plot(gps.mlt_apf9(i1),x,'.-'); grid on; datetick
	subplot(212); plot(gps.mlt_apf9(i1),y,'.-'); grid on; datetick
	disp('showing "good" GPS (from raw) in timeseries format');
	  disp('paused'); pause
	end

        elseif length(i1)==1,
	   % What to do if only one good GPS position (rare, but possible)?
	   % Try interpolating...
	   t0 = gps.mlt_apf9(i1);
	   lon0 = gps.lon(i1);
	   lat0 = gps.lat(i1);
	   if gotuni,
	     gd = find(isfinite(uni.got_gps) & isfinite(uni.u_sfc) & ...
		 isfinite(uni.v_sfc) & uni.got_gps & uni.flid==fid);
	     %gd = find(isfinite(uni.got_gps) & isfinite(uni.u_sfc) & isfinite(uni.v_sfc));
	     %gd = gd(find(uni.got_gps(gd) & uni.flid(gd)==fid));
	     if ~isempty(gd),
	       %u_sfc(ii) = sum(uni.u_sfc([ii-1 ii+1]))/2;
	       u_sfc(ii) = interp1(uni.utc_gps(gd),uni.u_sfc(gd),utc_gps(ii),'linear');
	       v_sfc(ii) = interp1(uni.utc_gps(gd),uni.v_sfc(gd),utc_gps(ii),'linear');
	       %v_sfc(ii) = sum(uni.v_sfc([ii-1 ii+1]))/2;
	       cx = [u_sfc(ii) 0];
	       cy = [v_sfc(ii) 0];
	     end
	   end
	end % if length(i1)>=2,


	% Now compute ubs,vbs from summing up accumulated velocities since
	% last surface interval, then purge variables

	if isfinite(lat_desc_vbs) % descent position has been saved
	  if isnan(utc_up(ii))
	    warning('Got GPS data but didn''t compute utc_up. Why?');
	  end
	
	    % ugps,vgps are from GPS descend-to-up displacement divided by time
	    % uef,vef are from time integrated Uacc,Vacc divided by the same time	

	    subsurface_time = (utc_up(ii) - utc_desc_vbs)*86400;
	    for i = 1:nmiss,
	      % remove surface drift time from integration
	      if argo_mode(upinds(i)),
		subsurface_time = subsurface_time - (utc_down(upinds(i)+1)-utc_up(upinds(i)))*86400;
	      else
		subsurface_time = subsurface_time - (utc_down(dninds(i))-utc_up(upinds(i)))*86400;
	      end
	    end
	    latR = (lat_up(ii) + lat_desc_vbs)/2;
	    Delx_gps = (lon_up(ii) - lon_desc_vbs)*lat2m*cos(latR*pi/180);
	    Dely_gps = (lat_up(ii) - lat_desc_vbs)*lat2m;

	    pmiss = find(isnan(Pacc));
	    if ~isempty(pmiss),
	      warning(['(Missing ' num2str(length(pmiss)) ...
		    ' pressure values)']);
	    end

	    % ugps,vgps are no longer used	    
	    ugps = Delx_gps/(utc_up(ii) - utc_desc_vbs)/86400;
	    vgps = Dely_gps/(utc_up(ii) - utc_desc_vbs)/86400;

	    Delx_sfc = 0; Dely_sfc = 0;
	    for i = 1:nmiss,
	      if argo_mode(upinds(i)),
		Delx_sfc = Delx_sfc + u_sfc(upinds(i))*(utc_down(upinds(i)+1)-utc_up(upinds(i)));
		Dely_sfc = Dely_sfc + v_sfc(upinds(i))*(utc_down(upinds(i)+1)-utc_up(upinds(i)));
	      else
		Delx_sfc = Delx_sfc + u_sfc(upinds(i))*(utc_down(dninds(i))-utc_up(upinds(i)));
		Dely_sfc = Dely_sfc + v_sfc(upinds(i))*(utc_down(dninds(i))-utc_up(upinds(i)));
	      end
	    end
	    
	    % So definition "vef = v - vbs" requires the removal of surface
	    % drift from vgps to get v (subsurface) comparable to vef
	    ugsub = (Delx_gps - Delx_sfc)/subsurface_time;
	    vgsub = (Dely_gps - Dely_sfc)/subsurface_time;

	    %	    u_gps(ii) = ugps;
	    %	    v_gps(ii) = vgps;
	    u_gps(ii) = ugsub;
	    v_gps(ii) = vgsub;

	    if ~ctdonly & ~isempty(Uacc)
	    % Modify Timacc to include down and up surface points (assuming
	    % no shear between surface and first measurement)
	    Timacc = [0; Timacc; utc_up(ii)-utc_desc_vbs];
	    Uacc = [Uacc(1); Uacc; Uacc(end)];
	    Vacc = [Vacc(1); Vacc; Vacc(end)];
	    Pacc = [0; Pacc; 0];
	    if length(Uacc)~=length(Timacc),
	      warning('Uacc and Timacc don''t match!');
	    end

	    [xv,Uacc] = fillnan(Timacc,Uacc);
	    [xv,Vacc] = fillnan(Timacc,Vacc);
	    
	    Delx_ef = cumtrapz(Timacc*86400,Uacc); % m
	    Dely_ef = cumtrapz(Timacc*86400,Vacc); % m
	
	    % Here make some plots showing the horizontal displacements and the
	    % pressure timeseries behavior and data locations, along with the
	    % surface drift...
	    
	    uef = Delx_ef(end)/subsurface_time;
	    vef = Dely_ef(end)/subsurface_time;
	
	    
	    %if ii==96
	    %  error('hit 96');
	    %end
	
	    
%	    ubs(ii) = ugps - uef;
%	    vbs(ii) = vgps - vef;
	    ubs(ii) = ugsub - uef;
	    vbs(ii) = vgsub - vef;

	    % and "corrected" Uacc,Vacc
	    Delx_efc = cumtrapz(Timacc*86400,Uacc+ubs(ii)); % m
	    Dely_efc = cumtrapz(Timacc*86400,Vacc+vbs(ii)); % m

	    % Now, fill in missing ubs,vbs for all previous profiles in vbs
	    % interval and u_sfc,v_sfc for previous missed profiles
	    if ~argo_mode(ii)
	      if nmiss>0,
		ubs([dninds(1)-2 dninds]) = ubs(ii);
		vbs([dninds(1)-2 dninds]) = vbs(ii);
		ubs(upinds) = ubs(ii);
		vbs(upinds) = vbs(ii);
	      else
		ubs(ii-1) = ubs(ii);
		vbs(ii-1) = vbs(ii);
	      end
	    else
	      if nmiss>0,
		ubs(upinds) = ubs(ii);
		vbs(upinds) = vbs(ii);
	      end
	    end
	    
	    if do_vbsplots
	      figure(1)
	      clf
	      subplot(321)
	      plot(Timacc,Pacc,'bx',...
		  [0 utc_up(ii)-utc_desc_vbs],[0 0],'ro');
	      axis ij; grid on; ylabel('pressure');
	      title('o: utc\_desc\_vbs, utc\_up(ii); x: Timacc');
	      subplot(323)
	      plot(Timacc,Delx_ef,'bx',...
		  Timacc,Delx_efc,'m.',...
		  [0 Timacc(end)],[0 Delx_gps],'ro--');
	      grid on; ylabel('x (m)');
	      subplot(325)
	      plot(Timacc,Dely_ef,'bx',...
		  Timacc,Dely_efc,'m.',...
		  [0 Timacc(end)],[0 Dely_gps],'ro--');
	      grid on; ylabel('y (m)');
	      xlabel('days');
	      subplot(122); plot(Delx_ef,Dely_ef,'.-',...
		  Delx_efc,Dely_efc,'m.',...
		  [0 Delx_gps],[0 Dely_gps],'ro--');
	      grid on; axis equal
	      title(['EM-APEX ' num2str(fltid) '; hpid:' num2str(hprof)...
	    '; ssid:' num2str(ssid(ii))]);
	      disp(['ugsub,vgsub = ' num2str(ugsub) ',' num2str(vgsub)]);
	      disp(['uef,vef = ' num2str(uef) ',' num2str(vef)]);
	      disp(['ubs,vbs = ' num2str(ubs(ii)) ',' num2str(vbs(ii))]);
	      disp(['usfc,vsfc = ' num2str(u_sfc(ii)) ',' num2str(v_sfc(ii))]);
	    
	      figure(2)
	      subplot(121)
	      plot(Uacc,Pacc,'r.-',Uacc+ubs(ii),Pacc,'m-');
	      %hold on; plot(u_sfc(ii),0,'k*',u_sfc(ii-2-nmiss*2),0,'cs');
	      %hold on; plot(u_sfc(ii),1,'k*',u_sfc(ii-2-nmiss*2),1,'cs');
	      hold on; plot(u_sfc(ii),1,'k^');
	      plot(Uacc(1)+ubs(ii),Pacc(1),'mv',Uacc(end)+ubs(ii),Pacc(end),'m^');
	      if vbsid(ii)>0
		plot(u_sfc(ii-2-nmiss*2),1,'kv');
	      end
	      if nmiss>=1
		plot(u_sfc(upinds),upinds*0,'go');
	      end
	      yl = get(gca,'ylim');
	      plot([0 0],yl,'k-', ubs(ii)*[1 1],yl,'k--');
	      plot(ugsub*[1 1],[0 max(Pacc)],'g-',uef*[1 1],[0 max(Pacc)],'g--');
	      %title('g-: ubar, g--: uefbar, k--: ubs, k*: usurf-aft, cs: usurf-bef');
	      title('g-: ubar, g--: uefbar, k--: ubs, k\^: usurf-aft, kv: usurf-bef');
	      grid on; axis ij; xlabel('r: Uacc; (m: "absolute")');
	      ylabel('pressure'); hold off
	      
	      subplot(122)
	      plot(Vacc,Pacc,'b.-',Vacc+vbs(ii),Pacc,'c-');
	      %hold on; plot(v_sfc(ii),0,'k*',v_sfc(ii-2-nmiss*2),0,'ms');
	      %hold on; plot(v_sfc(ii),1,'k*',v_sfc(ii-2-nmiss*2),1,'ms');
	      hold on; plot(v_sfc(ii),1,'k^');
	      plot(Vacc(1)+vbs(ii),Pacc(1),'cv',Vacc(end)+vbs(ii),Pacc(end),'c^');
	      if vbsid(ii)>0,
		plot(v_sfc(ii-2-nmiss*2),1,'kv');
	      end
	      if nmiss>=1
		plot(v_sfc(upinds),upinds*0,'go');
	      end
	      yl = get(gca,'ylim');
	      plot([0 0],yl,'k-', vbs(ii)*[1 1],yl,'k--');
	      plot(vgsub*[1 1],[0 max(Pacc)],'g-',vef*[1 1],[0 max(Pacc)],'g--');
%	      title('g-: vbar, g--: vefbar, k--: vbs, k*: vsurf-aft, ms: vsurf-bef');
	      title('g-: vbar, g--: vefbar, k--: vbs, k\^: vsurf-aft, kv: vsurf-bef');
	      hold off
	      grid on; axis ij; xlabel('b: Vacc (c: "absolute")');
	      disp('Showing accumulated pressure and position timeseries (v1)');
	      if pause_after
		disp('paused'); pause
	      end
	    end
	  end % if ~ctdonly
	    

	else % 	if isfinite(lat_desc)
	  % Profile has GPS data (surfacing) but not previous position
	  % (descent). For example, the first surfacing of a deployment (drop pos not
	  % known) or the first good fix after surface intervals with no GPS.
	  disp(['ubs,vbs not computed for float ' num2str(fid) ' profile ' num2str(hpid(ii))...
		' because no information on descent position.']);
	  % This looks like the place to fill in missing surface data since
	  % the last good GPS fix ("multi-surface module"). Calc up and down positions, then abs
	  % velocity params, by assuming a constant
	  % subsurface drift between the last known descent and the current
	  % surfacing plus surface displacements derived from
	  % intervening (interpolated) u_sfc,v_sfc and utc_up,utc_down values.
	  	
	  % Identify index of last known down position. lat_down should only be
	  % defined for a) descents from the surface or b) Argo mode ascents (in
	  % which case it represents the previous (missing) descent position.
	  if gotuni&isfield(uni,'got_gps'),
	    %idn = max(find(isfinite(lat_down) & uni.flid==fid & uni.hpid<hprof &...
	    %	[0 uni.got_gps(1:end-1)] ));
	    %idn = max(find(isfinite(lat_down) & uni.flid==fid & uni.hpid<=hprof));
	    idn = max(find(isfinite(lat_down) & flid==fid & hpid<=hprof));
	  else
	    idn = [];
	  end
	  
	end % if isfinite(lat_desc)
	
	% Re-initialize variables
	Uacc = [];
	Vacc = [];
	Timacc = [];
	Pacc = [];
	[utc_desc_vbs,lat_desc_vbs,lon_desc_vbs] = deal(NaN);

	if verbose
	  disp(['Finished vbarstar interval ' num2str(vbsnum)...
		' at profile ' num2str(hpid(ii))]);
	end

	
	vbsnum = vbsnum + 1; % vbs interval is an interval between good GPS fixes

	nmiss = 0; upinds = []; dninds = [];
	
	got_gps(ii) = 1;
      else % if ~isempty(i1),
	got_gps(ii) = 0;
      end % if ~isempty(i1),
      
      if isfinite(utc_up(ii)) & got_gps(ii)==0 & gotuni & isfield(uni,'got_gps')
	% If GPS didn't accomplish what we wanted, see whether gaps can be
	% filled in by interpolating from the previously saved unifiles
	% output. got_gps will remain 0 for this profile to indicate that 
	% variables (u_sfc, lat_up, lat_down, ubs) are all made up.

	% first step: fill in u_sfc and v_sfc for this (up) profile
	% simplest is to use linear interpolation, but maybe try something
	% fancier later...
	
	gd = find(isfinite(uni.got_gps) & isfinite(uni.u_sfc) & ...
	    isfinite(uni.v_sfc) & isfinite(uni.utc_up) & uni.got_gps & uni.flid==fid);
	%gd = find(isfinite(uni.got_gps) & isfinite(uni.u_sfc) & isfinite(uni.v_sfc) & isfinite(uni.utc_up));
	%gd = gd(find(uni.got_gps(gd) & uni.flid(gd)==fid));
	%if ~isempty(gd),
	if length(gd)>=2,
	u_sfc(ii) = interp1(uni.utc_up(gd),uni.u_sfc(gd),utc_up(ii),'linear');
	v_sfc(ii) = interp1(uni.utc_up(gd),uni.v_sfc(gd),utc_up(ii),'linear');
	end
	% next step: fill in lat_up,lon_up by assuming a constant subsurface
	% drift between the last known descent and the next known surfacing
	% plus surface displacements derived from intervening u_sfc,v_sfc and
	% utc_up,utc_down values. Actually, better to wait and do this
	% altogether at the next good GPS fix.

	% project Uacc,Vacc,Pacc,Timacc to the surface position, then drop to
	% zero velocity during the surface interval. later, when displacement
	% is calculated, we will add back in the surface drifts (which are
	% absolute, not EF, velocity and so are incompatible until we know
	% vbs)
	if ~ctdonly & nvel(ii)>nvmin,
	  Timacc = [Timacc; utc_up(ii)-utc_desc_vbs; utc_up(ii)-utc_desc_vbs+1/86400];
	  Uacc = [Uacc; Uacc(end); 0];
	  Vacc = [Vacc; Vacc(end); 0];
	  Pacc = [Pacc; 0; 0];
	  if gotuni & length(uni.utc_down)>=ii+1,
	    % add on end of surface interval using utc_down and shallowest
	    % velocity from next profile
	    Timacc = [Timacc; uni.utc_down(ii+1)-utc_desc_vbs; uni.utc_down(ii+1)-utc_desc_vbs+1/86400];
	    [minp,ish] = min(uni.Pef(:,ii+1));
	    Uacc = [Uacc; 0; uni.U(ish,ii+1)];
	    Vacc = [Vacc; 0; uni.V(ish,ii+1)];
	    Pacc = [Pacc; 0; 0];
	  end
	end
	% Note: the above "paragraph" adds points to the accumulated
	% variables that are probably redundant with the main surface
	% extrapolations done earlier. Is this a problem? As long as time
	% values used are identical, the duplicates shouldn't contribute to
	% the integral, although the routine might complain.
	
	% Just bump up count of missing surface GPS intervals for later use
	nmiss = nmiss + 1;
	upinds = [upinds ii];
	if gotuni & length(uni.argo_mode)>=ii+1 & isfinite(uni.argo_mode(ii+1)) & ~uni.argo_mode(ii+1),
	  dninds = [dninds ii+1];
	end
      end % if isfinite(utc_up(ii)) & got_gps(ii)==0 & gotuni & isfield(uni,'got_gps')

      
      % Continue to estimate the maximum number of depths needed
      maxvals = max(maxvals,npvals);
      if maxvals>Ndeps,
	  warning([num2str(maxvals) ' CTD values needed but only ' ...
		num2str(Ndeps) ' spaces allocated. Fix unifiles.']);
	end
      if ~ctdonly
	maxefvals = max(maxefvals,nefpvals);
	if maxefvals>Nefdeps,
	  warning([num2str(maxefvals) ' EF values needed but only ' ...
		num2str(Nefdeps) ' spaces allocated. Fix unifiles.']);
	end
      end
        
      if pause_each
	disp(['Paused at end of profile ' num2str(hpid(ii)) ' processing.']); pause
      end
	
      if surfaced(ii) & dnup(ii)==2,
	if verbose
	  disp(['Finished subsurface interval at profile ' num2str(hpid(ii))]);
	end

	[utc_desc,lat_desc,lon_desc] = deal(NaN);
      end
      ii = ii+1;
    end % if npvals<=2
  end %  if ~uponly | mod(aa(3),2)==0,
end %   for iD = 1:nD,

clear cx cy t0 lat0 lon0 ogps oi1

if 1
  if ctdonly
    disp(['Float ' num2str(fid) ' deployment ' num2str(deploy) ': ' ...
	  num2str(length(hpids)) ' half profiles, ' num2str(maxvals) ...
	  ' p levels max']);
  else
    disp(['Float ' num2str(fid) ' deployment ' num2str(deploy) ': ' ...
	  num2str(length(hpids)) ' half profiles, ' num2str(maxvals) ...
	  ' p levels max, ' num2str(maxefvals) ' Pef levels max.']);
  end
end


if do_plots_after
  iif = find(flid==fid&depl==deploy);
  dotsize = 5;
  %ax = [datenum(2005,7,30,20,0,0) datenum(2005,08,31,0,0,0) 0 150];
  %ax = [min(utc_gps(iif)) max(utc_gps(iif)) 0 560];
  %ax = [min(min(UTC(:,iif))) max(max(UTC(:,iif))) 0 560];
  ax = [min(min(UTC(:,iif))) max(max(UTC(:,iif))) 0 2000];

  % plot time-depth sections for U,V,T,S
  %subplot(411)
  subaxis2(4,1,1)
  scatter(squash(UTC(:,iif)),squash(P(:,iif)),dotsize,squash(T(:,iif)),'filled'); axis ij;
  hold on;
  axis(ax)
  datetick('x',6,'keeplimits');
  grid on
  %caxis([18.75 28]);
  %colormap jet; cb1 = colorbar;
  colormap jet(32); cb1 = colorbar;
  set(get(cb1,'ylabel'),'string','Temperature');
  title(fltid)
  ylabel('Depth (m)');
  %box on
  subaxis2(4,1,2)
  %scatter(ctd_mlt(:),P(:),dotsize,sigth(:),'filled'); axis ij;
  scatter(squash(UTC(:,iif)),squash(P(:,iif)),dotsize,squash(S(:,iif)),'filled'); axis ij;
  hold on;
  axis(ax)
  datetick('x',6,'keeplimits');
  grid on
  %caxis([36.4 36.9]);
  colormap jet; cb2 = colorbar;
  set(get(cb2,'ylabel'),'string','Salinity');
  ylabel('Depth (m)');
  box on
  subaxis2(4,1,3)
%  scatter(squash(UTCef(:,iif)),squash(Pef(:,iif)),dotsize,squash(U1(:,iif)),'filled'); axis ij;
  scatter(squash(UTCef(:,iif)),squash(Pef(:,iif)),dotsize,squash(U(:,iif)),'filled'); axis ij;
  hold on;
  %tdh = plot(mlt,pdep,'k.');
  axis(ax)
  datetick('x',6,'keeplimits');
  grid on
  %caxis([-1 1]);
  colormap jet; cb3 = colorbar;
  set(get(cb3,'ylabel'),'string','U');
  ylabel('Depth (m)');
  box on
  subaxis2(4,1,4)
%  scatter(squash(UTCef(:,iif)),squash(Pef(:,iif)),dotsize,squash(V1(:,iif)),'filled'); axis ij;
  scatter(squash(UTCef(:,iif)),squash(Pef(:,iif)),dotsize,squash(V(:,iif)),'filled'); axis ij;
  hold on;
  axis(ax)
  datetick('x',6,'keeplimits');
  grid on
  %caxis([-1 1]);
  colormap jet; cb4 = colorbar;
  set(get(cb4,'ylabel'),'string','V');
  %set(get(cb4,'ylabel'),'string','V2');
  ylabel('Depth (m)');
  xlabel('Date');
  %box on
  %pause
end

  if pause_after
    pause
  end

end % for fltidind=1:length(fltids),

if do_save
  if ctdonly
    if use_v6
      %outfile = [outdir 'allctd-v6.mat'];
    
      save(outfile,'P','T','S','UTC','utc_gps','hpid','ssid','lat_gps',...
	  'lon_gps','flid','depl','ind','utc_dep','utc_rec','-v6');
    else
      %outfile = [outdir 'allctd.mat'];
      
      save(outfile,'P','T','S','UTC','utc_gps','hpid','ssid','lat_gps',...
	  'lon_gps','flid','depl','ind','utc_dep','utc_rec',...
	  'maxp','botdep','u_sfc','v_sfc','u_gps','v_gps',...
	  'got_gps','apf9_terr');
    end
  else % if ctdonly
    if uponly
      %outfile = [outdir 'allprofs_up.mat'];
    else
      %outfile = [outdir 'allprofs.mat'];
    end
    save(outfile,'P','T','S','UTC','utc_gps','hpid','ssid','lat_gps',...
      'lon_gps','flid','depl','ind','utc_dep','utc_rec',...
      'maxp','botdep','vbsid','surfaced','argo_mode',...
      'U','V','Pef','UTCef','magvar',...
      'U1','V1','U2','V2','nvel',...
      'Verr1','Verr2','V1woW','V2woW','Vchan',...
      'utc_up','utc_down','ubs','vbs','u_sfc','v_sfc','u_gps','v_gps',...
      'got_gps','lat_up','lon_up','lat_down','lon_down',...
      'gps_dist','apf9_terr',...
      'Wp','Wr','Wf','ppos',...
      'Wp_ca','ppos_ca','P_ca','mlt_ca',...
      'ar','Nskip');
end

if exist('webdir')&exist(webdir)
  cmd = ['cp -f ' outfile ' ' webdir];
  unix(cmd)
end
end

% Notes (4/2/09): -- Integrated velocity for VBS computation assumes no shear
%   between shallowest measurement and surface (except for patching on
%   up-profile). Another option would be to use the surface drift for the
%   surface value. Need to make sure that ubs,vbs isn't added to surface
%   value.
% -- Accumulated velocity is always U1,V1 channel (rotated into Ug,Vg). Need
%   to add a check that this channel is better (or at least as good) on each
%   profile. [fixed, 3/10/10]
