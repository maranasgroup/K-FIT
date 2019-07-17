File1 = 'InputFile1.xlsx';
File2 = 'InputFile2.xlsx';
File3 = 'InputFile3.xlsx';

model=modcompile(File1, File2, File3);
res=kineticestimate(model,[]);
