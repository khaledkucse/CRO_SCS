package ABC;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;

public class Main 
{
	public static void main(String[] args) 
	{
		args = new String[6];
		//Input file destination
		args[0] = new String("Input/Random_Set/Small/4_5_10_7.txt");
		//Type: DNA or Protein
		args[1]="D";
		//No Of Bees
		args[2]="100";
		//Max Cycle Number
		args[3]="100";
		//Output file destination
		args[4]= new String("Output/Random_Set/Small/4_5_10_7_ABC_SCS_Result.txt");
		//Real Dataset || Random Dataset
		args[5] = "Random";
		
		if(args.length<4)
		{
			System.out.println("ERROR : Insufficient argument!!!");
			return;
		}
		
		String inputFile = args[0];	
		int status=0;
		Vector<Vector<Character>> input = new Vector<Vector<Character>>();
		String n;
		int num_strings= 0;
		int no_of_bees,max_cycle_number;
		double tStart= 0;
		double tEnd= 0;
		double Run_time=0;
		boolean valid =true;
		int size=0;
		Random r = new Random();
		try
		{
			Scanner sc =new Scanner(new File(inputFile));
			while(sc.hasNextLine())
			{
				n=sc.nextLine();
				Vector<Character> v = new Vector();
				for(int i=0;i<n.length();i++)
				{
					v.add(n.charAt(i));
				}
				input.add(v);
				
				num_strings++;
				size=v.size();
			}
			

		}
		
		catch(FileNotFoundException E)
		{
			System.out.println("ERROR : input file not found ! Please enter file's path and name as the first argument.");
			return;
		}
		//Start Parameter Check
		if(!(args[1].equalsIgnoreCase("D")||args[1].equalsIgnoreCase("P")))
		{
			System.out.println("ERROR : argument 2 must be \"D\" (for real DNA sequences), \"P\" (for real Protein sequences)!!");
			return;
		}
		
		
		try
		{
			no_of_bees= Integer.parseInt(args[2]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 2 is for Number of bees.It must be a integer value!!!");
			return;
		}
		try
		{
			max_cycle_number= Integer.parseInt(args[3]);
			//valid = false;
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 3 is for Maximum Cycle Number.It must be a Integer value!!!");
			return;
		}
		
		// end parameter check
		
		//output file name
		String OutputFile;
		if(args.length==6 && !(args[4].equals("")))
		{
			OutputFile = args[4];
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
			OutputFile = temp + "_ABC_SCS_Result.txt";
		}
		
		//Artificial Bee colony will run here
		Vector<Double> fitness = new Vector<Double>(no_of_bees);
		tStart= System.currentTimeMillis();
		Function func =new Function(input,no_of_bees,args[1]);
		//func.initializeSearching(num_strings,size,args[5],args[1]);
		if(valid==true)
		{
			func.frequencyCalculation();
			func.generateSCS(size);
			for(int i=0;i<max_cycle_number;i++)
			{
				
				fitness = func.fitnessCalculation();
				
				double avg = 0;
				for(int j=0;j<fitness.size();j++)
				{
					avg = avg+fitness.elementAt(j);
				}
				avg = avg/fitness.size();
				
				int rq=0;
				for(int j=0;j<fitness.size();j++)
				{
					if(fitness.elementAt(j)<avg)
					{
						rq=1;
						func.generateSCS(j,size);
					}
				}
				if(rq==0)
					break;
			}
		}
		
		
		tEnd= System.currentTimeMillis();	
		Run_time=(tEnd-tStart);
		
		//Dataset Name Retrival
		String[] token = inputFile.split("/");
		String temp=token[token.length-1];
		String[] hash= temp.split("t");
		temp=hash[0].substring(0,(hash[0].length()-1));
		
		//print the output
		String output = "";
		Vector<Character> v = func.showResult(num_strings,size,args[5],args[1]);
		output="Dataset: "+temp+" \n"+"SCS length= " + v.size() + "\nRun-time= "+ Run_time + " milliseconds" ;
		System.out.println(output);			
		output=output+ ";\nThe resulting SCS is: ";	
		
		for(int i=0;i<v.size();i++)
		{
			output=output+v.elementAt(i);
		}

		//print the resulting SCS to the output file
		int ret=Myprint(output,OutputFile);
		if(ret==1)
		{
			System.out.println("The resulting SCS was printed in : " + OutputFile);
		}
	}
	
	
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

}
