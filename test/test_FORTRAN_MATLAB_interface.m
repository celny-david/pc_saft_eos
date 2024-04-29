%% treat the .so as C library
clc
clear

% libpath = "/home/dixiw/DEV/CLionProjects/PC-SAFT/pc_saft/make/library_dynamic/";
libname = "pcsaft";

tic
if not(libisloaded(libname))
%     addpath(libpath)
    [notfound, warnings] = loadlibrary("libpcsaft.so", "libpcsaft.h",'alias',libname);
end
if not(libisloaded(libname))
    error("Library "+libname+" was not loaded");
end
libfunctions(libname,'-full')

console_file = "fort.6";

fclose(fopen(console_file,'wt')); % clear content of the file
t_init = toc;
% === call sequence ===
calllib('pcsaft','parse_input');

calllib(libname,'set_density',900.0)
calllib(libname,'set_temperature',350.0)
calllib(libname,'set_molar_ratio',[0.2,0.8])

p = libpointer('doublePtr',0.0);
pp = libpointer('doublePtr',0.0);
ppp = libpointer('doublePtr',0.0);

calllib(libname,'press_calc',p,pp,ppp)


p.Value
pp.Value
ppp.Value

% calllib('pcsaft','parse_input')

unloadlibrary(libname)
t_calc = toc;

disp("total time: "+num2str(t_calc))
disp("press calc: "+num2str(t_calc - t_init))
% === display the fortran console output ===
disp(fileread(console_file))

% %% treat the .so as c++ library
% 
% clibgen.buildInterface('libpcsaft.h')
