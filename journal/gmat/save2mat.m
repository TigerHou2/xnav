function save2mat(filepath,rData,vData,mu,i,j,k)
%SAVE2MAT saves pos/vel/mu data in a 3D cell structure from GMAT
%   Detailed explanation goes here

if ~isfile(filepath)
    temp = matfile(filepath,'Writable',true);
    temp.rArray = {};
    temp.vArray = {};
    temp.muArray = {};
end
variableInfo = who('-file', filepath);

if ismember('rArray', variableInfo), load(filepath,'rArray');
else, rArray = {}; end
if ismember('vArray', variableInfo), load(filepath,'vArray');
else, vArray = {}; end
if ismember('muArray',variableInfo), load(filepath,'muArray');
else, muArray = {};end

rArray{i,j,k}  = vclean(rData);
vArray{i,j,k}  = vclean(vData);
muArray{i,j,k} = mu;
save(filepath,'rArray','vArray','muArray');

end

