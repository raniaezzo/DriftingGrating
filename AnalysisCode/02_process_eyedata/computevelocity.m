%% this function is used to compute velocity
function v = computevelocity(d,samplingRateData)
         v = zeros(size(d));    

         %%calculate mean velocity for x,y seperately, and then sampliing
         % (window = 3)
         % multiplying by 0.5 to result in the average between:
         % <-- 2 and <-- 1 difference; and --> 1 and --> 2 difference
         % then multiply that by by a sampling rate /3 (window = 3,
         % because d(5:end,:) - d(2:end-3,:) is 3 samples offset and
         % d(4:end-1,:) - d(1:end-4,:) is also 3 samples offset.
         v(3:length(d)-2,:) = samplingRateData/3*(0.5*(d(5:end,:) - d(2:end-3,:) + d(4:end-1,:) - d(1:end-4,:)));

        %%calculate begin and end of segment, using window = 2

         % add points for the start of the velocity vector
         v(2,:) = samplingRateData/2*(d(3,:) - d(1,:));

         % add points for the end of the velocity vector
         v(length(d)-1,:) = samplingRateData/2*(d(end,:) - d(end-2,:));
end