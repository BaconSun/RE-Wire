function  RecoCurveCore(path0, path, outpath, group_name,options)
% exe_str=sprintf('%sCorresUI.exe %s %s %s %f %f %f %f %f %f %d %d %d',path0,path, outpath, group_name, options.select_thre,options.epiAreaThre,options.reguLamda,options.isoDScale,options.dire_thre,options.tangentThre,options.offsize,options.considerArea, options.comb);
exe_str=sprintf('%sreconstruct_3D_curves_node %s %s %s %f %f %f %f %f %f %d %d %d',path0,path, outpath, group_name, options.select_thre,options.epiAreaThre,options.reguLamda,options.isoDScale,options.dire_thre,options.tangentThre,options.offsize,options.considerArea, options.comb);
exe_str
system(exe_str);
end
