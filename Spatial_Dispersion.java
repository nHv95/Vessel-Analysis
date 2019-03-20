import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.plugin.*;
import ij.process.*;
import ij.measure.*;
import java.io.*;
import java.lang.*;
import java.awt.*;
import java.awt.event.*;
import java.lang.Math.*;
import java.awt.image.*;

/*
Given a set of isotropic, equidimensional binary stacks of squared images, this plugin calculates the number of 
cycles of rhombicubocthedric dilation necessary to fill 90%, 95% and 99% of each volume 
on the basis of the distribution of the black pixels. The plugin asks for a source folder containing
a set of subfolders (mandatory) containing image stacks with the same dimensions.
The procedure  is iterated over the nested folders for all the images present in every folder. 
The results are saved in a target location, choosed by the operator, as a number of .txt result files 
- one for each subfolder - named after the source subfolders.
Results are expressed as raw values (Hv) which are normalized, subfolder after subfolder,
according to the highest percentual volume observed in each subfolder to give nHv results. To 
normalize all subfolders in a single run, the image presenting the highest percentual volume 
should be present in all subfolders.

This plugin has been validated against simple images manually dilated and against complex images 
expanded by a macro ImageJ routine. This plugin requires at least ImageJ 1.43 and has been tested on 
version 1.44o. 

*/

public class Spatial_Dispersion implements PlugIn {

	public void run(String arg) {
 
// Setup name and version

		String thisplug = "Spatial_Dispersion";
		String versione = "1.0 del 30-04-2012";

// Choosing the folder to analyze and the target location for results

	   	String titolo1 = "Select the folder grouping all the folders in analysis";
	   	String titolo2 = "Where should I save the result files?";

		IJ.showMessage("Source", "Select the folder grouping all the folders in analysis");
		DirectoryChooser dc1 = new DirectoryChooser(titolo1);
		String dir1 = dc1.getDirectory();

		IJ.showMessage("Target", "Where should I save the result files?");
		DirectoryChooser dc2 = new DirectoryChooser(titolo2);
		String dir2 = dc2.getDirectory();

// Assessing starting time

		java.util.Date attuale = new java.util.Date();		String miaData = java.text.DateFormat.getDateTimeInstance(2, 2).format(attuale);
 
// Getting subfolders in input folder

		java.io.File sottodirs = new java.io.File(dir1);
		String[] lista = sottodirs.list();

// For every subfolder of the input folder...... 

		int sottocartella = 0;
		for (sottocartella=0; sottocartella<lista.length; sottocartella++) {
			String uc = lista[sottocartella].substring(0,1);
			if (uc.equals(".")) {
			}
			else {

// get the full path of the result file 

				String restit = dir2+lista[sottocartella]+".txt";

// getting the header of the file

				String identificativo = "Calculating normalized Hv indexes in nested folders.\nPlugin "+thisplug+" Version: "+versione+"\n\n";
				String firstline = "Analyzed folder of folders: "+dir1+"\nAt time: "+miaData+"\n\n";
				
// Creating a list of all subfolders

				String prosottodir1 = dir1+"/"+lista[sottocartella]+"/";
				java.io.File sottodir1 = new java.io.File(prosottodir1);
				String[] sublista = sottodir1.list();

// Defining occourring variables

				boolean Norm = false;
				double N = 0.00;
				double Ntmp = 0.00;
				int Cube = -1;
				int DimX = 0;
				int DimY = 0;
				int DimZ = 0;
				boolean primaImg = true; 
				int stackindex = 0;
				ImageStack stack = new ImageStack();
				double[][] Vini = new double [lista.length+1][sublista.length+1];

// For every stack image in the current subfolder

				for (stackindex=0; stackindex<sublista.length; stackindex++) {

					String suc = sublista[stackindex].substring(0,1);
					if (suc.equals(".")) {
					}
					else {

// open visible stacks one after the other

						Opener op =  new Opener();
						TiffDecoder prostack = new TiffDecoder(prosottodir1, sublista[stackindex]);
						try {
							FileInfo[] myinfo = prostack.getTiffInfo();
							ImagePlus imp = op.openTiffStack(myinfo);

// Check dimensions

							int[] dimensioni = imp.getDimensions();
							stack = imp.getStack();

							if (primaImg) {
								DimX = dimensioni[0];
								DimY = dimensioni[1];
								DimZ = dimensioni[3];
								primaImg = false;
							}
							else {
								if (dimensioni[0] != DimX) {
									IJ.error("Analysis Error", "Not all the stacks have the same width!");
								}
								if (dimensioni[1] != DimY) {
									IJ.error("Analysis Error", "Not all the stacks have the same heigth!");
								}
								if (dimensioni[3] != DimZ) {
									IJ.error("Analysis Error", "Not all the stacks have the same depth!");
								}
							}

// Get the maximal percentual volume 

							ByteProcessor bp = new ByteProcessor(dimensioni[0], dimensioni[1]);
							ImageProcessor ip = bp;
							double Vsgn = 0.00;
							double Vtot = DimX*DimY*DimZ;
							for (int s=1; s<dimensioni[3]+1; s++) {
								ip = stack.getProcessor(s);
								int[] histo = ip.getHistogram();
								Vsgn += histo[255];
							}
	
							Ntmp = 100*Vsgn/Vtot;
							if (Ntmp > N) {
								N=Ntmp;
								Cube = stackindex;
								Norm = true;
							}
							Vini[sottocartella][stackindex] = Ntmp;

// Stack closure

							imp.close();

						}
						catch (IOException e) {
							IJ.error("Scrivi file", e.getMessage());
							return;
						}
					}
				}

// Group information related to normalization to be saved on the result file 

				String norminfo = "Normalization Vol % : \t"+N+"\trelated to image stack:\t"+sublista[Cube]+"\n\n";
				String testata = "       Image:	Vol %: \tHv90%: \tHv95%: \tHv99%: \tnHv90%: \tnHv95%: \tnHv99%: \n";

// Begins the core elaboration to get the output result string  

				String corpo = "";

// For every visible subfolder in input folder begins the dilation cycle on all volumes, one at a time 

				for (stackindex=0; stackindex<sublista.length; stackindex++) {

					String suc = sublista[stackindex].substring(0,1);
					if (suc.equals(".")) {
					}
					else {

// Values of the Percentual Volume (VP) array are set to 100 and other variables are set to 0

						double[] VP = new double[510];
						for (int g=0; g<510; g++) {
							VP[g] = 100.00;
						}
						double VolPerc = 0.00000;
						double[] Vol = new double[510];
						int[] intVol = new int[510];
						boolean Hv1 = false;
						boolean Hv2 = false;
						boolean Hv3 = false;
						boolean Hv4 = false;
						double Hv90 = 0.00;
						double Hv95 = 0.00;
						double Hv99 = 0.00;
						double HvN = 0.00;
						double nHv90 = 0.00;
						double nHv95 = 0.00;
						double nHv99 = 0.00;
						int c = 0;
						int r = 0;
						int i = 0;
						int h= 0;
						int l = 0;
						int[] croce = new int[9];
						croce[0] = 0; 
						croce[1] = 1; 
						croce[2] = 0; 
						croce[3] = 1; 
						croce[4] = 1; 
						croce[5] = 1; 
						croce[6] = 0; 
						croce[7] = 1; 
						croce[8] = 0; 
						int[] quadrato = new int[9];
						for (int cnv=0; cnv<9; cnv++) {
							quadrato[cnv] = 1; 
						} 

// Open and duplicate the actual stack. The first is used as input and the second as output

						Opener op =  new Opener();
						TiffDecoder actualstack = new TiffDecoder(prosottodir1, sublista[stackindex]);
						try {
							FileInfo[] actualinfo = actualstack.getTiffInfo();
							ImagePlus impPari = op.openTiffStack(actualinfo);

							ImageStack pstack = impPari.getStack();
							ImageStack dstack = new ImageStack(DimX, DimY, impPari.getProcessor().getColorModel());

// Setup of image processors

							ByteProcessor wbp = new ByteProcessor(DimX, DimY);
							BinaryProcessor wbn = new BinaryProcessor(wbp);
							ImageProcessor wip = wbn;
							ByteProcessor tbp = new ByteProcessor(DimX, DimY);
							BinaryProcessor tbn = new BinaryProcessor(tbp);
							ImageProcessor tip = tbn;
							ByteProcessor bbp = new ByteProcessor(DimX, DimY);
							BinaryProcessor bbn = new BinaryProcessor(bbp);
							ImageProcessor bip = bbn;
							ByteProcessor dbp = new ByteProcessor(DimX, DimY);
							BinaryProcessor dbn = new BinaryProcessor(dbp);
							ImageProcessor dip = dbn;

							int dimension = (DimX)*(DimY);

							byte[] nopixels = new byte[dimension];
							byte[] propixels = new byte[dimension];
							byte[] inputpixels = new byte[dimension];
							byte[] workpixels = new byte[dimension];
							byte[] tworkpixels = new byte[dimension];
							byte[] bworkpixels = new byte[dimension];
							byte[] outputpixels = new byte[dimension];
							byte[] finpixels = new byte[dimension];
							byte[] swappixels = new byte[dimension];

// Filling the new void stack with slices and calculating the percentual volume of the input stack

							for (int sl=0; sl<DimZ; sl++) {
								dstack.addSlice("Vuoto", dip);
							}

							VP[0] = Vini[sottocartella][stackindex];

// Performing dilation only in the presence of some signal

							if (VP[0] > 0) {

// Expansion cycle

								while (VolPerc <99.5) {
									ij.IJ.showProgress(stackindex, sublista.length);
									r++;
									if (r >1) {
										r=0;
									}

// Cycle counter

									c++;

// For every slice i of the volume.....

									for (i=1; i<DimZ+1; i++) {
										h = i-1;
										l = i+1;

// Get images i-1 e i+1. Creates images if they do not exist

										if (h>0) {
											tip = pstack.getProcessor(h);
											inputpixels = (byte[]) pstack.getPixels(h);
											tworkpixels = inputpixels.clone();
										}
										else {
											tip.threshold(255);
											tworkpixels = nopixels.clone();
										}
										if (l<DimZ+1) {
											bip = pstack.getProcessor(l);
											inputpixels = (byte[]) pstack.getPixels(l);
											bworkpixels = inputpixels.clone();
										}	
										else {
											bip.threshold(255);
											bworkpixels = nopixels.clone();
										}

// Choice of the convolution scheme according to cycle number 

										if (r==0) {

// 3D cross expansion of image i-1

											tip.setPixels(tworkpixels);
											tip.threshold(0);

// 3D cross expansion of image i

											inputpixels = (byte[]) pstack.getPixels(i);
											workpixels = inputpixels.clone();

											dip.setPixels(workpixels);
											dip.convolve3x3(croce);
											dip.threshold(0);

// 3D cross expansion of image i+1

											bip.setPixels(bworkpixels);
											bip.threshold(0);

// Getting the final image

											dip.copyBits(tip,0,0,Blitter.OR);
											dip.copyBits(bip,0,0,Blitter.OR); 										}
										else {

// Square espansion of image i-1

											tip.setPixels(tworkpixels);
											tip.convolve3x3(quadrato);
											tip.threshold(0);

// Square espansion of image i

											inputpixels = (byte[]) pstack.getPixels(i);
											workpixels = inputpixels.clone();

											dip.setPixels(workpixels);
											dip.convolve3x3(quadrato);
											dip.threshold(0);

// Square espansion of image i+1

											bip.setPixels(bworkpixels);
											bip.convolve3x3(quadrato);
											bip.threshold(0);

// Getting the final image

											dip.copyBits(tip,0,0,Blitter.OR);
											dip.copyBits(bip,0,0,Blitter.OR);
										}

// Calculates the percentual volume of the actual image and save the result in an array

										int[] histo = dip.getHistogram();
										intVol[i] = histo[255];
										Vol[i] = 100*histo[255];

// The output stack is modified with the new images risulting from the actual voxel dilation 

										outputpixels = (byte[]) dip.getPixels();
										dstack.setPixels(outputpixels, i);
// End of the for cycle
									}	

// Stack swap!! Now PariDil becames Scarto, DispDil becames PariDil and Scarto is discarded

									for (int ls=1; ls<DimZ+1; ls++) {
										finpixels = (byte[]) dstack.getPixels(ls);
										swappixels = finpixels.clone();
										pstack.setPixels(swappixels, ls);
									}	

									for (int ls=1; ls<DimZ+1; ls++) {
										dstack.setPixels(nopixels, ls);
									}	

// Calculates the percentual volume of the actual image and save the result in an array

									double Vact = 0;
									for (int j=1; j<DimZ+1; j++) {
										Vact += intVol[j];
									}
									VolPerc = (Vact/(dimension*DimZ))*100;
									VP[c] = VolPerc;

// Calculates the percentual volume between 2 expansion cycles by linear interpolation

									if (r==0) {

										if (c<500) {
											VP[c-1] = (VP[c]+VP[c-2])/2;
										}
										else {
											VolPerc = 100;
										}
									}

// End of the while cycle

								} 

// End of instruction: if VP[0] > 0.

							}
							else {
								for (int g=1; g<510; g++) {
									VP[g] = 0.00;
								}
							}

// Getting raw Hv 3D indexes 

							for (int z=1; z<500; z++) {
								if (VP[z]>90) {
									while (Hv1 == false) {
										Hv90 = (z-1)+(90-VP[z-1])/(VP[z]-VP[z-1]);
										Hv1 = true;
									}
								}
							}
							for (int z=1; z<500; z++) {
								if (VP[z]>95) {
									while (Hv2 == false) {
										Hv95 = (z-1)+(95-VP[z-1])/(VP[z]-VP[z-1]);
										Hv2 = true;
									}
								}
							}
							for (int z=1; z<500; z++) {
								if (VP[z]>99) {
									while (Hv3 == false) {
										Hv99 = (z-1)+(99-VP[z-1])/(VP[z]-VP[z-1]);
										Hv3 = true;
									}
								}
							}
							if (Hv1 == false) {
								Hv90 = -1;
							}
							if (Hv2 == false) {
								Hv95 = -1;
							}
							if (Hv3 == false) {
								Hv99 = -1;
							}

// Getting normalized Hv 3D indexes if VP[0] > 0

							if (VP[0] == N) {
								HvN = 0;
							}
							if (VP[0] < N) {
								if (VP[0] > 0) {
									for (int z=1; z<500; z++) {
										if (VP[z] > N) {
											while (Hv4 == false) {
												HvN = (z-1)+(N-VP[z-1])/(VP[z]-VP[z-1]);
												Hv4 = true;
											}
										}
									}
								}
							}
							if (Hv4 == false) {
								HvN = 0;
							}

// Closing the stack

							impPari.close();

// End of the dilation cycles on a volume - let's start on a new one

						}
						catch (IOException e) {
							IJ.error("Scrivi file", e.getMessage());
							return;
						}

// Grouping result information on raw and normalized data 

						String startoutput = sublista[stackindex];
						float fVP = (float) VP[0];
						float fHvN = (float) HvN;
						float fHv90 = (float) Hv90;
						float fHv95 = (float) Hv95;
						float fHv99 = (float) Hv99;
						if (Norm) {
							nHv90 = Hv90-HvN;
							nHv95 = Hv95-HvN;
							nHv99 = Hv99-HvN;
							float fnHv90 = (float) nHv90;
							float fnHv95 = (float) nHv95;
							float fnHv99 = (float) nHv99;
							corpo = corpo.concat(startoutput+"\t"+java.lang.Float.toString(fVP)+"\t"+java.lang. Float.toString(fHv90)+"\t"+java.lang. Float.toString(fHv95)+"\t"+java.lang. Float.toString(fHv99)+"\t"+java.lang. Float.toString(fnHv90)+"\t"+java.lang. Float.toString(fnHv95)+"\t"+java.lang. Float.toString(fnHv99)+"\n");
						}
 						else {
							corpo = corpo.concat(startoutput+"\t"+java.lang. Float.toString(fVP)+"\t"+java.lang. Float.toString(fHv90)+"\t"+java.lang. Float.toString(fHv95)+"\t"+java.lang. Float.toString(fHv99)+"\t\n");
						}
 					}
				}

// Collect end of elaboration time for speed considerations 

				java.util.Date finale = new java.util.Date();				String datafinale = "\nFine elaborazione alle: " + java.text.DateFormat.getTimeInstance(2).format(finale);

// Opens the output file, writes the results and closes the file

				try {
					BufferedWriter risultati = new BufferedWriter(new FileWriter(restit));
					risultati.write(identificativo);
					risultati.write(firstline);
					risultati.write(norminfo);
					risultati.write(testata);
					risultati.write(corpo);
					risultati.write(datafinale);
					risultati.close();
				}
				catch (Exception e) {
					IJ.error("Scrivi file", e.getMessage());
					return;
				}

// End of the else instruction related to a subfolder  

			}

// End of operations on a subfolder - let's start witha new one

		}

		IJ.showMessage("End", "End of elaboration");
	}
}


