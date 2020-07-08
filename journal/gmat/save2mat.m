function save2mat(filepath,rData,vData,mu,i,j,k)
%SAVE2MAT saves pos/vel/mu data in a 3D cell structure from GMAT
%   Detailed explanation goes here

variableInfo = who('-file', filepath);

if ismember('rArray', variableInfo), load(filepath,'rArray');
else, rArray = {}; end
if ismember('vArray', variableInfo), load(filepath,'vArray');
else, vArray = {}; end
if ismember('muArray',variableInfo), load(filepath,'muArray');
else, muArray = {};end

rArray{i,j,k}  = rData;
vArray{i,j,k}  = vData;
muArray{i,j,k} = mu;
save(filepath,'rArray','vArray','muArray');

end

