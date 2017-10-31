function [output tasks] = setup_NL(model)
    
        tasks = [];
        
        options.name='bike';
        options.select_thre=1.5; %The distance threshold for prunning curves in view3
        options.dire_thre= 0.8; %The tangential threshold for prunning curves in view3
        options.epiAreaThre= 8.0; %The neighboring scale of the degenerate case
        options.reguLamda= 50000; %The weight of the regularization term in the fitting for the degenerate case, generally it's set as 50000
        options.isoDScale= 1.0; %The ratio of the distance of the isolated curve and its nearest curve to the curve's length, generally it's set as 1
        options.offsize= 10; %The minimum length of a 3D curve candidate, generally it's set as 10
        options.tangentThre = 1.0; %A parameter in detecting the tangent point of the image curve, generally it's set as 1.0
        options.considerArea=20; %The neighboring scale in detecting the tangent point of the image curve, generally it's set as 20
        options.comb=1; %The index of the combinations of three views
        tasks = [tasks; options];

        options.name='human';
        options.select_thre=3.8;
        options.dire_thre= 0.8;
        options.epiAreaThre= 10.0;
        options.reguLamda= 50000;
        options.isoDScale= 0.5;
        options.offsize= 15;
        options.tangentThre = 1;
        options.considerArea=30;
        options.comb=5;
        tasks = [tasks; options];

        options.name='bird';
        options.select_thre=1.0;
        options.dire_thre= 0.5;
        options.epiAreaThre= 6.0;
        options.reguLamda= 80000;
        options.isoDScale= 0.5;
        options.offsize= 2;
        options.tangentThre = 0.2;
        options.considerArea=20;
        options.comb=4;
        tasks = [tasks; options];

        options.name='elephant';
        options.select_thre=1.0;
        options.dire_thre= 0.8;
        options.epiAreaThre= 7.5;
        options.reguLamda= 50000;
        options.isoDScale= 0.5;
        options.offsize= 12;
        options.tangentThre = 1.0;
        options.considerArea=20;
        options.comb=6;
        tasks = [tasks; options];

        options.name='fatcat';
        options.select_thre=1.5;
        options.dire_thre= 0.7;
        options.epiAreaThre= 9.0;
        options.reguLamda= 50000;
        options.isoDScale= 0.8;
        options.offsize= 5;
        options.tangentThre = 1;
        options.considerArea=30;
        options.comb=2;
        tasks = [tasks; options];

        options.name='horse';
        options.select_thre=2.5;
        options.dire_thre= 0.9;
        options.epiAreaThre= 5.0;
        options.reguLamda= 50000;
        options.isoDScale= 1.3;
        options.offsize= 10;
        options.tangentThre = 1.0;
        options.considerArea=10;
        options.comb=1;
        tasks = [tasks; options];
       
        options.name='flower';
        options.select_thre=1.5;
        options.dire_thre= 0.72;
        options.epiAreaThre= 7.0;
        options.reguLamda= 50000;
        options.isoDScale= 2.0;
        options.offsize= 5;
        options.tangentThre = 1;
        options.considerArea=30;
        options.comb=3;
        tasks = [tasks; options];

        options.name='cart';
        options.select_thre=1.5;
        options.dire_thre= 0.5;
        options.epiAreaThre= 9.0;
        options.reguLamda= 50000;
        options.isoDScale= 0.5;
        options.offsize= 5;
        options.tangentThre = 1.0;
        options.considerArea=20;
        options.comb=6;
        tasks = [tasks; options];

        options.name='turtle';
        options.select_thre=2.3;
        options.dire_thre= 0.85;
        options.epiAreaThre= 4.0;
        options.reguLamda= 50000;
        options.isoDScale= 1.0;
        options.offsize= 8;
        options.tangentThre = 1.0;
        options.considerArea=40;
        options.comb=2;
        tasks = [tasks; options];
        
        output = '';
        for i = 1: length(tasks)
            if strcmp(tasks(i).name, model)
                output = tasks(i);
            end
        end
end

