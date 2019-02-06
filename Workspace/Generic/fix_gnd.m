
id_gen = RandSampleTagGeneratorX(1:3);
ids = id_gen.generate(3);
timestamp = datetime();

N = size(GroundTruth,2);
velocity = zeros(1,2);
GroundTruthTracks = GroundTruthTrackX.empty(0,3);
for i = 1:N
    GroundTruth_i = cell2mat(GroundTruth(:,i));
    n = size(GroundTruth_i,2);
    clear GroundTruthList;
    if i<N
        velocity = cell2mat(GroundTruth(:,i+1))-GroundTruth_i;
    end
    for j = 1:n
        gndVector = [GroundTruth_i(1,j);
                         velocity(1);
                         GroundTruth_i(2,j);
                         velocity(2);];
        if i == 1
            GroundTruthTracks(j) = GroundTruthTrackX(GroundTruthStateX(gndVector,... 
                                                                       timestamp),...
                                                     ids(j));
        else
            GroundTruthTracks(j).Trajectory(end+1) = GroundTruthStateX(gndVector,... 
                                                                       timestamp, ...
                                                                       ids(j));
        end
    end
    timestamp = timestamp + duration(0,0,1);
    %GroundTruthListSequence{i} = GroundTruthList;
end

gndTruthTrack = GroundTruthTrackX(GroundTruthStateX, TagX(12434));