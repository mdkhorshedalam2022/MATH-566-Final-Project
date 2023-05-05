function Nodes(type)
    global M N nodetype
    
    for ii=N+2:(M-2)*N+N-1
        nodetype(ii)=0;%'Core';
    end
    
    for ii=1
        nodetype(ii)=1; %'bottom left node'
    end
    for ii=N
        nodetype(ii)=2; %'bottom right node'
    end
    
    for ii=(M-1)*N+1
        nodetype(ii)=3;%'top left node'
    end
        
    for ii=M*N
        nodetype(ii)=4; %'top right node'
    end
    
    
    
    for ii=2:N-1
        nodetype(ii)=5;%'bottom';
    end
    
    for ii= (M-1)*N+2:M*N-1
        nodetype(ii)=6;%'top'; 
    end
    
    
    
    for ii=N+1:N:(M-2)*N+1
        nodetype(ii)=7;%'left side'; 
    end
    
    for ii=2*N:N:(M-2)*N+N
        nodetype(ii)=8;%'right side';
    end
    
    if strcmp(type,'A')
        % Case A Horziontal Gradient 
        % A one-to-three Combiantion 
        %Midpoint left inlet
        for ii=(floor(M/2)-1)*N+1:N:(floor(M/2))*N+1
           nodetype(ii)=9; %inlet;
        end
        
        %Midpoint right outlet
        for ii=(floor(M/2)-1)*N+N:N:(floor(M/2))*N+N
           nodetype(ii)=10;%oulet;
        end
        %Endpoints right outlet
        for ii=[2*N, 3*N,(M-2)*N+N, (M-3)*N+N]
           nodetype(ii)=10;%oulet;
        end
    
    end
end