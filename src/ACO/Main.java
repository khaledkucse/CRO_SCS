package ACO;

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

		args = new String[11];
		//Input file destination
		args[0] = new String("Input/Random_Set/Small/4_5_10_7.txt");
		//Type: DNA or Protein
		args[1]="D";
		//No Of Colony
		args[2]="2";
		//No Of Ant in each Colony 
		args[3]="16";
		//q0
		args[4]="0.6";
		//p
		args[5]="0.5";
		//Alpha
		args[6]="9";
		//Beta
		args[7]="9";
		//Gamma
		args[8]="1000";
		//Output file destination
		args[9]= new String("Output/Random_Set/Small/4_5_10_7_ACO_SCS_Result.txt");
		//Real Dataset || Random Dataset
		args[10] = "Real";
		
		if(args.length<9)
		{
			System.out.println("ERROR : Insufficient argument!!!");
			return;
		}
		
		String inputFile = args[0];	
		int status=0;
		Vector<Vector<Character>> input = new Vector<Vector<Character>>();
		Vector<Vector<Character>> supersequence = new Vector<Vector<Character>>();
		Vector<Double> pheromone =new Vector<Double>();
		String n;
		int num_strings= 0;
		int no_of_colony,no_of_ant,alpha,beta,gamma;
		double q0,p;
		int VT = 0;
		double tStart= 0;
		double tEnd= 0;
		Vector<Character> o =new Vector<Character>();
		double Run_time=0;
		boolean valid =true;
		int size=0;
		Random r = new Random();
		
		try
		{
			Scanner sc =new Scanner(new File(inputFile));
			System.out.println("Reading Data from File...");
			while(sc.hasNextLine())
			{
				n=sc.nextLine();
				Vector<Character> v = new Vector();
				//Vector<Integer> b = new Vector();
				for(int i=0;i<n.length();i++)
				{
					v.add(n.charAt(i));
				}
				input.add(v);
				
				num_strings++;
				size=v.size();
			}
			
			System.out.println("Executing ACO algorithm...");
			tStart= System.currentTimeMillis();
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
			no_of_colony= Integer.parseInt(args[2]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 2 is for population size.It must be a integer value!!!");
			return;
		}
		try
		{
			no_of_ant= Integer.parseInt(args[3]);
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 3 is for no_of_ant.It must be a Integer value!!!");
			return;
		}
		try
		{
			q0=Double.parseDouble(args[4]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 4 is for q0.It must be a double value!!!");
			return;
		}
		try
		{
			p= Double.parseDouble(args[5]);	
			//valid=false;
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 5 is for Initial KE.It must be a double value!!!");
			return;
		}
		try
		{
			alpha= Integer.parseInt(args[6]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 6 is for KE Loss Rate.It must be a Double value!!!");
			return;
		}
		try
		{
			beta= Integer.parseInt(args[7]);	
		}
		catch(NumberFormatException E)
		{
			System.out.println("ERROR : Argument 7 is for Decompotion reaction threshold.It must be a integer value!!!");
			return;
		}
		try
		{
			gamma= Integer.parseInt(args[8]);	
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
			if(num_strings == 5000 && size == 1000)
			{
				Function func1 = new Function(input,1);
				func1.initializeSearching(num_strings,size,"Random", "D");
				o = func1.showResult(); 
				valid = func1.getTrail(3);
			}
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
			OutputFile = temp + "_ACO_SCS_Result.txt";
		}
		
		//ACO Run
		
		Function func =new Function(input,testcolony());
		Function colony2 = new Function(input,testcolony());
		func.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
		colony2.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
		//func.initializeSearching(num_strings,size,args[10],args[1]);
		
		if(valid == true)
		{
			Vector<Vector<Double>> prevpheromone = initializepheromone(num_strings,size);
			Vector<Character> z = new Vector<Character>();
			for(int iterate =0; iterate<50;iterate++)
			{
				for(int j=1;j<=8;j++)
				{
					Vector<Character> v =new Vector<Character>();
					v = func.generateSCS(); 
					z = colony2.generateSCS();
					supersequence.add(v);
					if(z != null)
						supersequence.add(v);
					double x = func.updateTrails(j);
					pheromone.add(x);
					x = colony2.updateTrails(j); 
					if(x !=0)
						pheromone.add(x);
					
					func.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
					colony2.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
				}
				int w = elitAnt(supersequence);
				func.setSupersequence(supersequence.elementAt(w));
				colony2.setSupersequence(supersequence.elementAt(w));
				for(int j=9;j<=16;j++)
				{
					Vector<Character> v =new Vector<Character>();
					v = func.generateSCS(); 
					z = colony2.generateSCS();
					supersequence.add(v);
					if((z != null))
						supersequence.add(v);
					double x = func.updateTrails(j);
					pheromone.add(x);
					x = colony2.updateTrails(j); 
					if(x !=0)
						pheromone.add(x); 
					func.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
					colony2.initialize(input, no_of_colony, no_of_ant, q0, p, alpha, beta, gamma);
				}
				Vector<Vector<Double>> prespheromone = func.getPheromone();
				
				for(int i=0;i<prespheromone.size();i++)
				{
					Vector<Double> v = prevpheromone.elementAt(i);
					Vector<Double> s = prespheromone.elementAt(i);
					for(int j=0;j<s.size();j++)
					{
						double temp = s.elementAt(j)-v.elementAt(j);
						temp = (1-p) * s.elementAt(j)+gamma * temp;
						s.set(j, temp);
					}
					prespheromone.set(i, s);
				}
				
				func.setPheromone(prespheromone);
				prevpheromone = prespheromone;
				
			}
			
			//for(int i=0;i<80;i++)
			//{
				//System.out.println(supersequence.elementAt(i).size());
			//}
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
		Vector<Character> v = new Vector<Character>();
		if(valid == true)
			v = func.showResult(supersequence);
		else
			v=o;
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
	private static int elitAnt(Vector<Vector<Character>> supersequence)
	{
		int min =0;
		for(int i=0;i<supersequence.size();i++)
		{
			if(supersequence.elementAt(i).size() < supersequence.elementAt(min).size())
				min =i;
		}
		return min;
		
	}
	static int flag =0;
	private static int testcolony()
	{
		if(flag==0)
		{
			flag=1;
			return 1;
		}
		return 2;
	}
	private static Vector<Vector<Double>> initializepheromone(int num_string, int length)
	{
		Vector<Vector<Double>> result = new Vector<Vector<Double>>(num_string);
		for(int i=0;i<num_string;i++)
		{
			Vector<Double> v = new Vector<Double>();
			for(int j=0;j<length;j++)
			{
				v.add(1.00);
			}
			result.add(v);
		}
		
		return result;
		
	}

}
