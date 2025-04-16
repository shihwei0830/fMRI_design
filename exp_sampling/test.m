opt = {'CollectOutput',true};
hdr = {};
out = {};
[fid,msg] = fopen('data.txt','rt');
assert(fid>=3,msg) % ensure the file opened correctly.
fgetl(fid); % read and ignore the very first line.
while ~feof(fid)
    hdr{end+1} = fgetl(fid);
    out(end+1) = textscan(fid,'%f%f',opt{:});
end
fclose(fid);