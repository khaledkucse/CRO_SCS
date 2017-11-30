package DRandIBS_SCS;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;



public class FindSCS 
{
	public static void main(String[] args) throws IOException 
	{
		//args for IBS_SCS = InputFile string_type IBS_SCS/DR beam_size num_best OutputFile(optional)
		//args for DR =      InputFile string_type IBS_SCS/DR OutputFile(optional)
		args=new String[6];
		args[0]=new String("Input/Random_Set/Small/4_5_10_7.txt");
		args[1]="D";
		args[2]="IBS_SCS";
		args[3]="100";
		args[4]="7";
		args[5]=new String("Output/Random_Set/Small/4_5_10_7_IBS_SCS_Result.txt");
		if(!(args.length==3 || args.length==4 || args.length==5 ||args.length==6 ))
		{
			System.out.println("ERROR : invalid number of arguments! Please read \"readme.txt\"");
			return;
		}

		String inputFile = args[0];		
		int num_strings= 0;
		try
		{
			LineNumberReader reader  = new LineNumberReader(new FileReader(inputFile));				
			while (( reader.readLine()) != null)
			{
				num_strings++;
			}		
			reader.close();

		}
		catch(FileNotFoundException E)
		{
			System.out.println("ERROR : input file not found ! Please enter file's path and name as the first argument.");
			return;
		}

		if(!(args[1].equalsIgnoreCase("D")||args[1].equalsIgnoreCase("P")||args[1].equalsIgnoreCase("O")))
		{
			System.out.println("ERROR : argument 2 must be \"D\" (for real DNA sequences), \"P\" (for real Protein sequences), or \"O\" (otherwise)!");
			return;
		}
		if(!(args[2].equalsIgnoreCase("DR") || args[2].equalsIgnoreCase("IBS_SCS")))
		{
			System.out.println("ERROR : argument 3 must be either DR or IBS_SCS!");
			return;
		}
		if(args[2].equalsIgnoreCase("DR")&& !(args.length==3 || args.length==4))
		{
			System.out.println("ERROR : too many arguments for DR. Please read \"readme.txt\"");
			return;
		}
		if(args[2].equalsIgnoreCase("IBS_SCS")&& !(args.length==5 || args.length==6))
		{
			System.out.println("ERROR : too few arguments for IBS_SCS. Please read \"readme.txt\"");
			return;
		}


		String dataFormat = args[1];
		String whichAlgorithm = args[2];
		String OutputFile;	
		long tStart,tEnd;
		double Run_time;

		if(whichAlgorithm.equalsIgnoreCase("IBS_SCS"))// IBS_SCS
		{
			//get beam size
			String beam_size_s = args[3];
			int beam_size;			
			try
			{
				beam_size= Integer.parseInt(beam_size_s);	
			}
			catch(NumberFormatException E)
			{
				System.out.println("ERROR : argument 4 for IBS_SCS is beam size and must be an integer!");
				return;
			}
			if(beam_size<=0)
			{
				System.out.println("ERROR : argument 4 for IBS_SCS is beam size and must be a positive integer!");
				return;
			}

			//get number of best dominators
			String num_best_s = args[4];
			int num_best;
			try
			{
				num_best= Integer.parseInt(num_best_s);
			}
			catch(NumberFormatException E)
			{
				System.out.println("ERROR : argument 5 for IBS_SCS is the number of best dominators and must be an integer!");
				return;
			}
			if(num_best<0)
			{
				System.out.println("ERROR : argument 5 for IBS_SCS is the number of best dominators and must be a nonnegative integer!");
				return;
			}
				

			//determine output filename
			if(args.length==6)
			{
				OutputFile = args[5];
			}
			else
			{
				String temp;
				if(inputFile.lastIndexOf('.')!= -1)
					temp = inputFile.substring( 0 , inputFile.lastIndexOf('.') );
				else
				{
					temp = inputFile;
				}
				OutputFile = temp + "_IBS_SCS_Result.txt";
			}

			//run IBS_SCS
			IBS_SCS objIBS_SCS= null;
			tStart= System.currentTimeMillis();	   	    
			objIBS_SCS=new IBS_SCS(num_strings,inputFile,dataFormat);
			objIBS_SCS.mainMethod(beam_size,num_best);
			tEnd= System.currentTimeMillis();	
			Run_time=(tEnd-tStart);	
			
			//Dataset Name Retrival
			String[] token = inputFile.split("/");
			String temp=token[token.length-1];
			String[] hash= temp.split("t");
			temp=hash[0].substring(0,(hash[0].length()-1));
			
			//print the output
			String output = "";
			output="Dataset: "+temp+" \n"+"SCS length= " + objIBS_SCS.retLen + "\nRun-time= "+ Run_time + " milliseconds" ;
			System.out.println(output);			
			output=output+ ";\nThe resulting SCS is: ";
			for(int i=0; i<objIBS_SCS.retLen; i++)
			{
				output = output + (objIBS_SCS.Alphabet.charAt( objIBS_SCS.Global_finalSol.X[i] ) ) ;
			}		

			//print the resulting SCS to the output file
			int ret=Myprint(output,OutputFile);
			if(ret==1)
			{
				System.out.println("The resulting SCS was printed in : " + OutputFile);
			}

		}

		if(whichAlgorithm.equalsIgnoreCase("DR"))//DR
		{
			//determine output filename
			if(args.length==4)
			{
				OutputFile = args[3];
			}
			else
			{
				String temp;
				if(inputFile.lastIndexOf('.')!= -1)
					temp = inputFile.substring(0 , inputFile.lastIndexOf('.') );
				else
				{
					temp = inputFile;
				}
				OutputFile = temp + "_DR_Result.txt";
			}

			//run DR
			DR objDR= null;
			tStart= System.currentTimeMillis();	   	    
			objDR=new DR(num_strings,inputFile,dataFormat);
			objDR.mainMethod();
			tEnd= System.currentTimeMillis();	
			Run_time=(tEnd-tStart);	

			//print the output
			String output = "";
			output="SCS length= " + objDR.retLen + "\nrun-time= "+ Run_time + " milliseconds" ;
			System.out.println(output);			
			output=output+ "; the resulting SCS is: ";
			for(int i=0; i<objDR.retLen; i++)
			{
				output = output + (objDR.Alphabet.charAt( objDR.Global_finalSol.X[i] ) ) ;
			}		

			//print the resulting SCS to the output file
			int ret=Myprint(output,OutputFile);
			if(ret==1)
			{
				System.out.println();
				System.out.println("The resulting SCS was printed in : " + OutputFile);
				System.out.println();
			}

		}

	}//main


	public static int Myprint(String st, String Outputfile)
	{
		PrintWriter pr;	
		try
		{
			pr = new PrintWriter(new FileWriter(Outputfile, true));
			pr.println(st);
			pr.println();
			pr.close();

		}
		catch(FileNotFoundException e)
		{
			System.out.println("ERROR : output file not found ! Please enter file's path and name as the last argument (optional).");
			return -1;
		}	 
		catch(IOException e)
		{
			System.out.println("ERROR in printing the output to the file.");
			return -1;
		}	 
		return 1;
	}

}