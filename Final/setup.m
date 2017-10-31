function [output tasks options] = setup(model)
    
    tasks = [];

    task.model     = 'human';
    task.m         = 1;   % if m == 1, use TSP, otherwise, go for mTSP
    task.lambda    = [0.1, 1.0, 0.03];  % lamda-qo, mu-qo, mu-tsp
    task.tol       = [8.e-7, 2.e+6]; % For B-spline fitting
    task.defloat   = [60, 60, 0.004, 0.005]; % For removing isolated curves
    task.desupport = [0.004, 0.003, 0.2, 0.2, cos(pi/3), cos(pi/3)]; % For removing isolated curves
    task.cleantsp  = false;
    tasks = [tasks; task];

    task.model     = 'horse';
    task.m         = 1;
    task.lambda    = [0.2, 0.8, 0.03]; 
    task.tol       = [7.e-7, 1.e+5];
    task.defloat   = [50, 50, 0.04, 0.05];
    task.desupport = [0.03, 0.035, 0.05, 0.06, 0.5, 0.5];
    task.cleantsp  = false;
    tasks = [tasks; task];
 
    task.model     = 'flower';
    task.m         = 1;
    task.lambda    = [0.1, 1.0, 0.06]; 
    task.tol       = [7.e-7, 7.e+5];
    task.defloat   = [80, 80, 0.02, 0.015];
    task.desupport = [0.06, 0.01, 0.03, 0.07, 0.5, 0.0];
    task.cleantsp  = false;
    tasks = [tasks; task];
    
    
    task.model     = 'bird';
    task.m         = 1;
    task.lambda    = [0.1, 1.0, 0.05]; 
    task.tol       = [7.e-7, 1.e+5];
    task.defloat   = [8, 8, 0.2, 0.8];
    task.desupport = [0.015, 0.05, 0.3, 0.142, 0.0, 0.5];
    task.cleantsp  = false;
    tasks = [tasks; task];   
   

    task.model     = 'elephant';
    task.m         = 1;
    task.lambda    = [0.1, 1.0, 0.039]; 
    task.tol       = [7.e-7, 2.e+5];
    task.defloat   = [40, 40, 0.05, 0.05];
    task.desupport = [0.02, 0.04, 0.38, 0.2, 0.6, 0.5];
    task.cleantsp  = true;
    tasks = [tasks; task];    


    task.model     = 'bike';
    task.m         = 1;
    task.lambda    = [0.001, 0.05, 0.1]; 
    task.tol       = [7.e-7, 3.e+5];
    task.defloat   = [40, 80, 0.1, 0.06];
    task.desupport = [0.078, 0.07, 0.2, 0.2, 0.5, 0.5];
    task.cleantsp  = false;
    tasks = [tasks; task];

    
    task.model     = 'fatcat';
    task.m         = 1;
    task.lambda    = [0.01, 0.1, 0.015];  
    task.tol       = [7.e-7, 1.e+6];
    task.defloat   = [60, 80, 0.008, 0.01];
    task.desupport = [0.005, 0.01, 0.045, 0.03, 0.5, 0.5];
    task.cleantsp  = false;
    tasks = [tasks; task];
    

    task.model     = 'cart';
    task.m         = 2;
    task.lambda    = [0.01, 1, 0.01];  
    task.tol       = [7.e-7, 15.e+4];
    task.defloat   = [60, 60, 0.02, 0.03];
    task.desupport = [0.021, 0.022, 0.08, 0.08, 0.6, 0.6];
    task.cleantsp  = false;
    tasks = [tasks; task];
    

    task.model     = 'turtle';
    task.m         = 3;
    task.lambda    = [0.5, 0.1, 0.002];  
    task.tol       = [7.e-7, 1.e+4];
    task.defloat   = [100, 80, 0.026, 0.038];
    task.desupport = [0.01, 0.022, 0.1, 0.12, 0.0, 0.77];
    task.cleantsp  = false;
    tasks = [tasks; task];

    

    
    %% read-out
    output = '';
    for i = 1: length(tasks)
        if strcmp(tasks(i).model, model)
            output = tasks(i);
            output.suffix = '';
            if output.m > 1
                output.type = 'm';
            else
                output.type = 's';
            end
        end
    end
    
    options = setup_NL(model);
end

