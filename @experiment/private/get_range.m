function rl = get_range(indname)

switch(indname)
    case {'hma', 'hme'}, rl = [-0.03 0.03];
    case {'ae1','ae3','ae2','aef','aee','aep','pae'}, rl = [-30 30];    
    case {'aa1','aa3','aa2','aaf','aae','aap'}, rl = [0 50];
    case 'jra', rl = [0 3];
    case {'jer','jef','jei'}, rl = [0 4];
    case 'pax', rl = [-5 5]; 
    case 'pay', rl = [-5 5]; 
    case 'sim', rl = [0 2];
    case 'ddu',rl = [0 1];
    case 'adu', rl = [0 0.7];
    case 'dur', rl = [0 1.5];
    case 'lat', rl = [-0.03 0.03];
    case 'caf',rl=[-1 1];
    case 'fer',rl=[0 0.05];
    case 'rti',rl=[0 1];
    case 'pax', rl = [-5 5];
    case 'pay', rl = [-5 5];
    case 'sim', rl = [0 1.5];
    case 'asp',rl = [0. 0.5];
    case 'psp', rl = [0.2 0.7];
        
end