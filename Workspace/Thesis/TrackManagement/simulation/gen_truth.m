function truth= gen_truth(model)

%variables
truth.K= 200;                   %length of data/number of scans
truth.X= cell(truth.K,1);             %ground truth for states of targets  
truth.N= zeros(truth.K,1);            %ground truth for number of targets
truth.L= cell(truth.K,1);             %ground truth for labels of targets (k,i)
truth.track_list= cell(truth.K,1);    %absolute index target identities (plotting)
truth.total_tracks= 0;          %total number of appearing tracks

%target initial states and birth/death times
nbirths= 10;
wturn = 2*pi/180;

xstart(:,1)  = [ 1000+3.8676; -2.2; 1500-11.7457; -9.2; wturn/2 ];              tbirth(1)  = randi([1,truth.K-55]);     tdeath(1)  = randi([tbirth(1)+50,truth.K+1]);
xstart(:,2)  = [ -250-5.8857;  2; 1000+11.4102; -20; -wturn/8 ];                 tbirth(2)  = randi([1,truth.K-55]);     tdeath(2)  = randi([tbirth(2)+50,truth.K+1]);
xstart(:,3)  = [ -2000; 4; 250+6.7993; 1; -wturn/3 ];                         tbirth(3)  = randi([1,truth.K-55]);     tdeath(3)  = randi([tbirth(3)+50,truth.K+1]);
xstart(:,4)  = [ -1500; 10; -500; 20; 0 ];                                       tbirth(4)  = randi([1,truth.K-55]);     tdeath(4)  = randi([tbirth(4)+50,truth.K+1]);
xstart(:,5)  = [ 500-3.8676; -5; 0-11.0747; 11; wturn/5 ];                       tbirth(5)  = randi([1,truth.K-55]);     tdeath(5)  = randi([tbirth(5)+50,truth.K+1]);
xstart(:,6)  = [ -250+7.3806; 12; 1000-6.7993; -2; wturn/3 ];                 tbirth(6)  = randi([1,truth.K-55]);     tdeath(6)  = randi([tbirth(6)+50,truth.K+1]);
xstart(:,7)  = [ 1000; -10; 1500; -10; wturn/4 ];                                 tbirth(7)  = randi([1,truth.K-55]);     tdeath(7)  = randi([tbirth(7)+50,truth.K+1]);
xstart(:,8)  = [ 0; 10; 750; 0; -wturn/6 ];                                    tbirth(8)  = randi([1,truth.K-55]);     tdeath(8)  = randi([tbirth(8)+50,truth.K+1]);
xstart(:,9)  = [ 1000; -10; 1500; 10; wturn/9 ];                                 tbirth(9)  = randi([1,truth.K-55]);     tdeath(9)  = randi([tbirth(9)+50,truth.K+1]);
xstart(:,10)  = [ 250; -10; 750; 0; wturn/2 ];                                 tbirth(10)  = randi([1,truth.K-55]);     tdeath(10)  = randi([tbirth(10)+50,truth.K+1]);

%generate the tracks
truth.Logs = [];
Log.xV = NaN(2,truth.K);
for targetnum=1:nbirths
    truth.Logs{targetnum} = Log;
end
GroundTruthStateSequence = cell(1,truth.K);
GroundTruthTracks = GroundTruthTrackX.empty();
timestamp_init = datetime;
for targetnum=1:nbirths
    targetstate = xstart(:,targetnum);
    k = tbirth(targetnum);
    timestamp = timestamp_init + seconds(k);
    state = GroundTruthStateX(targetstate(1:4,1), timestamp); 
    truth.Logs{targetnum}.xV(:,k) = [targetstate(1,1);targetstate(3,1)];
    truth.X{k}= [truth.X{k} targetstate];
    truth.track_list{k} = [truth.track_list{k} targetnum];
    truth.N(k) = truth.N(k) + 1;
    track = GroundTruthTrackX(state);
    GroundTruthStateSequence{k}(end+1) = state;
    for k=tbirth(targetnum)+1:min(tdeath(targetnum),truth.K)
        timestamp = timestamp_init + seconds(k);
        targetstate = gen_newstate_fn(model,targetstate,'noiseless');
        truth.Logs{targetnum}.xV(:,k) = [targetstate(1,1);targetstate(3,1)];
        truth.X{k}= [truth.X{k} targetstate];
        truth.track_list{k} = [truth.track_list{k} targetnum];
        truth.N(k) = truth.N(k) + 1;
        state = GroundTruthStateX(targetstate(1:4,1), timestamp);
        track.Trajectory(end+1) = copy(state);
        GroundTruthStateSequence{k}(end+1) = copy(state);
    end
    GroundTruthTracks(end+1) = copy(track);
end
truth.GroundTruthTracks = GroundTruthTracks;
truth.GroundTruthStateSequence = GroundTruthStateSequence(~cellfun('isempty',GroundTruthStateSequence));
truth.total_tracks= nbirths;

