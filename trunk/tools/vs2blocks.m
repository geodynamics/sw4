%
% vs2blocks
%
%  First use readProfile(fil) to get the vs, dh and freq...
%
%        vs2blocks( vs, dh, freq )
%
%        Input: vs   - Average vs value in z
%               dh   - Coarse dh for grid (maps to z=zmax vs)
%               freq - Frequency to resolve

function vs2blocks(vs, dh, freq)
    
    % Constant
    ppw=10

    % Compute the minimum grid spacing required to resolve this source
    % Note:  This will vary from source type to source type, so for now
    % we're just going to apply this simple forumla for them all...
    nz = length(vs);
    
    h = zeros(size(vs));
    for i=1:nz
        h(i) = vs(i)/(ppw*freq);
    end
    
    % Original region spans the whole grid
    z1 = 0.0;
    
    % We take one pass through backwards and compute the block indices,
    % from coarsest to finest.  Then we take one more pass and enforce
    % constaints of the block sizes (i.e., they need to be at least
    % 3 coarse grid points wide...)
    blockEnds = zeros(size(vs));
      
    % We are going to compute blocks based on this coarse value
    rr = 1;
    blockEndCount = 0;
    for i=nz:-1:1
       % Walk backwards from slowest region to fastest
       if h(i) <= dh/rr
          z2 = z1;
          z1 = i*dh;
          blockEnds(i) = 1;
          disp(['End at: ' num2str(i*dh)])
          blockEndCount = blockEndCount + 1;
          rr = rr*2;
       end
    end   
    disp(['Maximum refinement ratio: ' num2str(rr)])
    disp(['Number of blocks pre constraints: ' num2str(blockEndCount)])
    
    % Enforce size constraints...
    % First block starts at 1
    blockStartIndex = 1;
    blockEndCount = 0;
    for i=2:nz
        if (blockEnds(i) == 1) && ((i-blockStartIndex) < 3)
            % need to be at least 3 coarse grid points..., so
            % push it out
            blockEnds(i) = 0;
            blockEnds(i+1) = 1;
        elseif (blockEnds(i) == 1)
            disp(['End at: ' num2str(i*dh)])
            blockEndCount = blockEndCount +1;
            blockStartIndex = i;
        end
    end
    disp(['Number of blocks post constraints: ' num2str(blockEndCount)])
    clf;
    plot(h);
    hold;
    % Now write out blocks with constraints...
    z1 = 0.0;
    for i=2:nz
        if blockEnds(i) == 1
           z2 = i*dh;
           plot(i, h(i), 'dr', 'LineWidth',2, 'MarkerSize', 8, 'MarkerFaceColor', [1 0 0]);
           disp(['refinement rr=' num2str(rr) ' z1=' num2str(z1) ' z2=' num2str(z2)])  
           z1 = z2;
           rr = rr/2;
        end
    end
    hold;
    
return