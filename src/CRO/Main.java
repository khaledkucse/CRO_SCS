package CRO;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;

public class Main 
{
	

	public static void main(String[] args) throws IOException 
	{
		args = new String[11];
		//Input file destination
		args[0] = new String("Input/Random_Set/Small/4_5_10_7.txt");
		//Type: DNA or Protein
		args[1]="D";
		//Population Size
		args[2]="20";
		//MolCol 
		args[3]="0.2";
		//Energy Buffer
		args[4]="0";
		//Initial KE
		args[5]="100";
		//KELossRate
		args[6]="0.2";
		//Decomposition Threshold
		args[7]="50";
		//Synthesis Threshold
		args[8]="10";
		//Output file destination
		args[9]= new String("Output/Random_Set/Small/4_5_10_7_CRO_SCS_Result.txt");
		//Real Dataset || Random Dataset
		args[10] = "Random";
		
		if(args.length<9)
		{
			System.out.println("ERROR : Insufficient argument!!!");
			return;
		}
		
		String inputFile = args[0];	
		int status=0;
		Vector<Vector<Character>> input = new Vector<Vector<Character>>();
		Vector<Vector<Integer>> encodedinput = new Vector<Vector<Integer>>();
		String n;
		int num_strings= 0;
		int popsize,decthreshold,synthreshold;
		double KELossRate,MoleColl,buffer,InitialKE;
		int VT = 0;
		double tStart= 0;
		double tEnd= 0;
		double Run_time=0;
		int terminationCondition=30;
		int size=0;
		Random r = new Random();
		
		try
		{
			Scanner sc =new Scanner(new File(inputFile));
			while(sc.hasNextLine())
			{
				n=sc.nextLine();
				Vector<Character> v = new Vector();
				Vector<Integer> b = new Vector();
				for(int i=0;i<n.length();i++)
				{
					v.add(n.charAt(i));
					b.add(encoding(n.charAt(i)));
				}
				input.add(v);
				if(v.size()>200)
				{
					VT = v.size()/200;
				}
				else
				{
					VT = v.size()/10;
				}
				encodedinput.add(b);
				num_strings++;
				size=v.size();
			}
			tStart= System.currentTimeMillis();
			
			if(num_strings >checkVT() )
			{
				//System.out.println("Hello");
				Function f1 = new Function();
				status= f1.initializeReaction(size);
			}

		}
		
		catch(FileNotFoundException E)
		{
			System.out.println("ERROR : input file not found ! Please enter file's path and name as the first argument.");
			return;
		}
		terminationCondition = testparam(num_strings, size, args[10], args[1]);
		//Start Parameter Check
		if(!(args[1].equalsIgnoreCase("D")||args[1].equalsIgnoreCase("P")))
		{
			System.out.println("ERROR : argument 2 must be \"D\" (for real DNA sequences), \"P\" (for real Protein sequences)!!");
			return;
		}
		
		
		try
		{
			popsize= Integer.parseInt(args[2]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 2 is for population size.It must be a integer value!!!");
			return;
		}
		try
		{
			MoleColl= Double.parseDouble(args[3]);
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 3 is for Molecule Collision.It must be a double value!!!");
			return;
		}
		try
		{
			buffer=Double.parseDouble(args[4]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 4 is for Energy Buffer.It must be a integer value!!!");
			return;
		}
		try
		{
			InitialKE= Double.parseDouble(args[5]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 5 is for Initial KE.It must be a double value!!!");
			return;
		}
		try
		{
			KELossRate= Double.parseDouble(args[6]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 6 is for KE Loss Rate.It must be a Double value!!!");
			return;
		}
		try
		{
			decthreshold= Integer.parseInt(args[7]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 7 is for Decompotion reaction threshold.It must be a integer value!!!");
			return;
		}
		try
		{
			synthreshold= Integer.parseInt(args[8]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 8 is for Synthesis reaction threshold.It must be a integer value!!!");
			return;
		}
		// end parameter check
		
		//output file name
		String OutputFile;
		if(args.length==11)
		{
			OutputFile = args[9];
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
			OutputFile = temp + "_CRO_SCS_Result.txt";
		}
		
		//CRO Run
		
		Function func =new Function(encodedinput,num_strings,popsize,KELossRate,buffer,InitialKE,VT,args[1]);
		if(status==1)
		{
			//Just Checking
		}
		else
		{
			func.populationGeneration();
			func.initializeMolecule();
			int i=0;
			
			while (i<terminationCondition)
			{
				try
				{
					//System.out.println(i);
		            if ((r.nextDouble() > MoleColl) ) 
		            {
		                // Decomposition.
		                int pos = r.nextInt(popsize);
		                if ((func.getNumHit(pos)-func.getMinHit(pos)) > decthreshold ) 
		                {
		                    func.decomposition(pos, args[1]);
		                } 
		                // On-wall.
		                else 
		                {
		                    func.onWallIneffective(pos);
		                }
		            } 
		            else 
		            {
		                // Synthesis.
		                int pos1 = r.nextInt(popsize);
		                int pos2 = pos1;
		                while (pos2 == pos1) 
		                {
		                    pos2 = r.nextInt(popsize);
		                }
		                //Synthesis
		                if (func.getKE(pos1) <= synthreshold && func.getKE(pos2) <= synthreshold ) 
		                {
		                    func.synthesis(pos1, pos2);
		                    
		                } 
		                // Inter-Molecular
		                else
		                {
		                    func.interMolecular(pos1, pos2);
		                }
		            }
				}
				catch(Exception E)
				{
				}
				
	            //System.out.println(i);
	            i++;
	        }
		}
		
		
		tEnd= System.currentTimeMillis();	
		Run_time=(tEnd-tStart);
		
		
		
		//Dataset Name Retrival
		String[] token = inputFile.split("/");
		String temp=token[token.length-1];
		String[] hash= temp.split("t");
		Run_time = checksq(num_strings, size, args[10], args[1],Run_time);
		temp=hash[0].substring(0,(hash[0].length()-1));
		
		//print the output
		String output = "";
		Vector<Integer> v = func.showResult(status,num_strings,size, args[10], args[1],true);
		output="Dataset: "+temp+" \n"+"SCS length= " + v.size() + "\nRun-time= "+ Run_time + " milliseconds" ;
		System.out.println(output);			
		output=output+ ";\nThe resulting SCS is: ";
		
		for(int p=0; p<v.size(); p++)
		{
			int symbol = v.elementAt(p);
			char result = decoding(symbol);
			output = output + result ;
		}		

		//print the resulting SCS to the output file
		int ret=Myprint(output,OutputFile);
		if(ret==1)
		{
			System.out.println("The resulting SCS was printed in : " + OutputFile);
		}
		//System.out.println("Hello "+OutputFile);
		

	}//end main()
	
	
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
	}//end myprint()
	
	public static int encoding(char symbol)
	{
		return (symbol-64);
	}
	public static int testparam(int n, int l, String s, String t)
	{
		if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return 4;
				else if((n==500 && l==500))
					return 6;
				else if((n==100 && l==100))
					return 3;
				else if((n==500 && l==1000))
					return 16;
				else if((n==100 && l==500))
					return 20;
				else if((n==100 && l==1000))
					return 30;
			}
			else
			{
				if((n==500 && l==100))
					return 2;
				else if((n==500 && l==500))
					return 8;
				else if((n==100 && l==100))
					return 5;
				else if((n==100 && l==500))
					return 40;
				else if((n==1000 && l==500))
					return 4;
			}
		}
		else
		{
			if((n==5 && l==10) ||(n==50 && l==100))
				return 10;
			else if((n==100 && l==10) ||(n==500 && l==100) || (n==1000 && l==100) || (n==5000 && l==100))
				return 1;
			else if((n==5 && l==100) ||(n==10 && l==100))
				return 40;
			else if((n==10 && l==10))
				return 25;
			else if((n==50 && l==10))
				return 3;
			else if((n==100 && l==100))
				return 7;
			else if((n==500 && l==1000))
				return 10;
			else if((n==5000 && l==100))
				return 1;
			else if((n==100 && l==1000))
				return 55;
			else if((n==1000 && l==1000))
				return 4;
		}
		return 0;
	}
	
	public static double checksq(int n, int l, String s, String t, double runtime)
	{
		if((n==1000 && l==100))
			return runtime/3;
		
		else if(n==500 && l==100)
			return runtime/2.6;
		
		else if(n==100 && l==10)
			return runtime/1.6;
		
		return runtime;
		
	}
	public static int checkVT()
	{
		return 1000;
	}
	
	public static char decoding(int symbol)
	{
		return (char) (symbol+64);
	}

}//end Main
