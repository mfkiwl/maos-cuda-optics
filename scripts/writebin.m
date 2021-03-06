%%Guide lines: 
%Each block of data can optionally be preceeded with a header
%If the data is cell array, the whole cell can have a header, and each entry in
%the cell can also have a header. %If they all have header, when the file is
%read out, the header of the whole cell will be appeneded to the end of the
%list of the headers for each cell. The same applies when trying to save the file.


function writebin(data,varargin)
    if nargin==3
        header=varargin{1};
        fn=varargin{2};
    elseif nargin==2
        fn=varargin{1};
        header='';
    else
        error('Usage: writebin(data,header,filename) or writebin(data, filename)');
    end
    do_gzip=0;%default don't do compress
              %if fn ends with .bin, we don't do gzip
    if length(fn)>3 &&strcmp(fn(end-3:end),'.bin')
        do_gzip=0;
    end
    %if fn ends with .gz, strip it off.
    if length(fn)>2 && strcmp(fn(end-2:end),'.gz')
        fn=fn(1:end-3);
        do_gzip=1;
    end
    %if fn doesn't end with .bin, add .bin
    if ~(length(fn)>3 &&strcmp(fn(end-3:end),'.bin'))
        do_gzip=0;
        fn=[fn '.bin'];
    end
    fid=fopen(fn,'wb');
    writebin_do(data,header,fid);
    if do_gzip
        system(sprintf('gzip -f %s',fn));
    end
    
    fclose(fid);
function writeheader_do(header, fid)
     M_HEADER=25856;
    if ~ischar(header)
        error('header must be char type');
    end
    fwrite(fid,uint32(M_HEADER),'uint32');
    fwrite(fid,uint64(length(header)+1),'uint64');
    fwrite(fid,header,'char');
    fwrite(fid,0,'char');
    fwrite(fid,uint64(length(header)+1),'uint64');
    fwrite(fid,uint32(M_HEADER),'uint32');
    
function writebin_do(data,header,fid)
    M_CSP64=25600;
    M_SP64=25601;
    M_CSP32=25606;
    M_SP32=25607;
    M_DBL=25602;
    M_INT64=25603;
    M_CMP=25604;
    M_INT32=25605;
    MCC_ANY=25633;
    MAT_SP=65281;
    MAT_CSP=65282;
    M_HEADER=25856;
    if length(size(data))>2
        error('We do not handle more than 2 dimensions');
    end
    if iscell(data)
        if ~isempty(header) 
            if ~iscell(header) %header is not cell.
                writeheader_do(header, fid);
            else
                if numel(header)==numel(data)+1
                    writeheader_do(header{numel(data)+1},fid);
                end
            end
        end
        fwrite(fid,MCC_ANY,'uint32');
        fwrite(fid,size(data,1),'uint64');
        fwrite(fid,size(data,2),'uint64');
        if ~isempty(header) && iscell(header)
            for ii=1:numel(data)
                writebin_do(data{ii},header{ii},fid);
            end
        else
            for ii=1:numel(data)
                writebin_do(data{ii},[],fid);
            end
        end
    else
        %write the header.
        if ~isempty(header)
            writeheader_do(header, fid);
        end
        if issparse(data)
            nz=nzmax(data);
            [Ir,Ic,P]=find(data);
            if isreal(P)
                complex=0;
            else
                complex=1;
                P2=zeros(2,length(P));
                P2(1,:)=real(P);
                P2(2,:)=imag(P);
                P=P2(:);
            end
            if nz<intmax('uint32') %use 32bit to save
                use_32=1;
                if complex
                    magic=M_CSP32;
                else
                    magic=M_SP32;
                end
            else
                use_32=0;
                if complex
                    magic=M_CSP64;
                else
                    magic=M_SP64;
                end
            end
            nx=size(data,1);
            ny=size(data,2);
            fwrite(fid,magic,'uint32');
            fwrite(fid,nx,'uint64');
            fwrite(fid,ny,'uint64');
            if nz>0
                fwrite(fid,nz,'uint64');
                Jc=zeros(ny+1,1);
                pos=1;
                for icol=1:ny
                    Jc(icol)=pos;
                    while pos<=nz && Ic(pos)==icol
                        pos=pos+1;
                    end
                end
                Jc(ny+1)=pos;
                %first convert to 0 index
                Ir=Ir-1;
                Jc=Jc-1;
                if use_32
                    Ir=uint32(Ir);
                    Jc=uint32(Jc);
                    fwrite(fid,Jc,'uint32');
                    fwrite(fid,Ir,'uint32');
                else
                    Ir=uint64(Ir);
                    Jc=uint64(Jc);
                    fwrite(fid,Jc,'uint64');
                    fwrite(fid,Ir,'uint64');
                end
                fwrite(fid,P,'double');
            end
        elseif isreal(data)
            fwrite(fid,M_DBL,'uint32');
            fwrite(fid,size(data,1),'uint64');
            fwrite(fid,size(data,2),'uint64');
            fwrite(fid,data,'double');
        else %complex
            fwrite(fid,M_CMP,'uint32');
            fwrite(fid,size(data,1),'uint64');
            fwrite(fid,size(data,2),'uint64');
            tmp=zeros(2,numel(data));
            tmp(1,:)=real(data(:));
            tmp(2,:)=imag(data(:));
            fwrite(fid,tmp,'double');
            clear tmp;
        end 
    end