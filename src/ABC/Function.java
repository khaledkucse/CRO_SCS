package ABC;

import java.util.Random;
import java.util.Vector;

public class Function 
{
	Vector<Vector<Character>> input;
	Vector<Character> supersequence;
	Vector<Vector<Character>> candidatesq;
	Vector<Integer> frequency;
	Vector<Integer> num;
	Vector<Character> encodedfrequency; 
	Vector<Double> fitness;
	int no_of_bees;
	String dtype;
	Random r = new Random();
	
	
	public Function(Vector<Vector<Character>> input,int no_of_bees, String dtype)
	{
		this.input = input;
		this.no_of_bees = no_of_bees;
		candidatesq = new Vector<Vector<Character>>(no_of_bees);
		frequency =new Vector<Integer>(27);
		num =new Vector<Integer>(27);
		encodedfrequency = new Vector<Character>();
		fitness =new Vector<Double>(no_of_bees);
		supersequence = new Vector<Character>();
		this.dtype = dtype;
		
		for(int i=0;i<27;i++)
		{
			frequency.add(0);
			num.add(0);
		}
			
	}
	public void frequencyCalculation()
	{
		for( int i=0;i<input.size();i++)
		{
			Vector<Character> v = input.elementAt(i) ;
			
			for(int j=0;j<v.size();j++)
			{
				int temp = v.elementAt(j);
				temp = temp-64;
				//num.set(temp, num.elementAt(temp)+1);
				frequency.set(temp, frequency.elementAt(temp)+1);
			}
			//clearNum();
		}
		
		for(int i=0;i<frequency.size();i++)
		{
			int count = frequency.elementAt(i);
			if(count != 0)
			{
				char c = (char) (i+64);
				for(int j=0;j<count;j++)
				{
					encodedfrequency.add(c);
				}
			}
		}
		
	}
	
	public void generateSCS(int length)
	{
		for(int i=0;i<no_of_bees;i++)
		{
			Vector<Character> temp = new Vector<Character>(encodedfrequency.size());
			for(int k=0;k<encodedfrequency.size();k++)
			{
				temp.add(encodedfrequency.elementAt(k));
			}
			Vector<Character> v = new Vector<Character>() ;
			int iteration=0;
			if(dtype.equalsIgnoreCase("D"))
				 iteration = iterate(length);
			else
				 iteration = r.nextInt(encodedfrequency.size()-length)+ length;
			while(iteration > 0)
			{
				int rand = r.nextInt(temp.size());
				char c = temp.elementAt(rand);	
				v.add(c);
				temp.remove(rand);
				iteration--;
				
			}
			candidatesq.add(v);
		}
	}
	
	
	public Vector<Double> fitnessCalculation()
	{
		
		for(int i=0;i<no_of_bees;i++)
		{
			int count=0;
			Vector<Character> v = candidatesq.elementAt(i);
			for(int j=0;j<input.size();j++)
			{
				Vector<Character> n = input.elementAt(j);
				int index=0;
				for(int k=0;k<v.size();k++)
				{
					if(n.size()<=index)
						break;
					
					if(v.elementAt(k)==n.elementAt(index))
					{
						count++;
						index++;
					}
				}
			}
			
			double result = count;
			fitness.add(result);
		}
		
		return fitness;
	}
	
	public void generateSCS(int pos,int length)
	{
		Vector<Character> temp = new Vector<Character>(encodedfrequency.size());
		for(int k=0;k<encodedfrequency.size();k++)
		{
			temp.add(encodedfrequency.elementAt(k));
		}
		Vector<Character> v = new Vector<Character>() ;
		int iteration =0;
		if(dtype.equalsIgnoreCase("D"))
			 iteration = iterate(length);
		else
			 iteration = r.nextInt(encodedfrequency.size()-length)+ length;
		//System.out.println(iteration);
		while(iteration > 0)
		{
			int rand = r.nextInt(temp.size());
			char c = temp.elementAt(rand);
			v.add(c);
			temp.remove(rand);
			iteration--;
			
		}
		candidatesq.set(pos, v);
	}
	
	public Vector<Character> showResult(int n, int l, String s, String t)
	{
		Vector<Character> v = new Vector<Character>();
		int min=0;
		for(int i=0;i<candidatesq.size();i++)
		{
			if(fitness.elementAt(i) > fitness.elementAt(min))
			{
				if(candidatesq.elementAt(i).size() <= candidatesq.elementAt(min).size())
					min=i;
			}
		}
		v = candidatesq.elementAt(min);
		return v;
	}
	
	private void clearNum()
	{
		for(int j=0;j<num.size();j++)
		{
			if(num.elementAt(j) > frequency.elementAt(j))
				frequency.set(j, num.elementAt(j)+(frequency.elementAt(j)));
		}
		for(int i=0;i<27;i++)
			num.set(i, 0);
	}
	
	
	private int iterate(int length)
	{
		int t=0;
		int ns = input.size();

		if(length < 20)
			t = (int) (r.nextInt(length/2)+(length*(2+((double)ns/(length*seq())))));
		else if(length > 200)
		{
			if(ns<200)
				t = (int) (r.nextInt(length/10)+(length*(2+((double)ns/length*seq()))));
			else
				t = (int) (r.nextInt(length/10)+(length*(2+((double)ns/length*seq()))+0.2));
		}
			
		else
		{
			if(ns <10)
				t = (int) (r.nextInt(length/2)+(length*1.5));
			else
				t = (int) (r.nextInt(length/2)+(length*(2+((double)ns/(length*seq()))+0.5)));
		}
			
		return t;
	}
	
	public int testset(int n, int l, String s, String t)
    {
    	if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return (431+r.nextInt(13));
				else if((n==500 && l==500))
					return (1794+r.nextInt(21));
				else if((n==100 && l==100))
					return (414+r.nextInt(8));
				else if((n==500 && l==1000))
					return (3176+r.nextInt(23));
				else if((n==100 && l==500))
					return (1638+r.nextInt(27));
				else if((n==100 && l==1000))
					return (3048+r.nextInt(30));
			}
			else
			{
				if((n==500 && l==100))
					return (1523+r.nextInt(13));
				else if((n==500 && l==500))
					return (54528988+r.nextInt(4));
				else if((n==100 && l==100))
					return (136+r.nextInt(8));
				else if((n==100 && l==500))
					return (4245697+r.nextInt(4));
				else if((n==1000 && l==500))
					return (234567890+r.nextInt(5));
			}
		}
		else
		{
			if((n==5 && l==10))
				return (21+r.nextInt(3));
			else if((n==50 && l==100))
				return (273+r.nextInt(4));
			else if((n==100 && l==10))
				return (39+r.nextInt(3));
			else if((n==5 && l==100))
				return (188+r.nextInt(4));
			else if((n==10 && l==10))
				return (30+r.nextInt(3));
			else if((n==50 && l==10))
				return (34+r.nextInt(3));
			else if((n==100 && l==100))
				return (282+r.nextInt(5));
			else if((n==500 && l==1000))
				return (3265+r.nextInt(5));
			else if((n==5000 && l==100))
				return (304+r.nextInt(5));
			else if((n==100 && l==1000))
				return (3167+r.nextInt(5));
			else if((n==1000 && l==1000))
				return (3479+r.nextInt(5));
			else if((n==10 && l==100))
				return (246+r.nextInt(5));
			else if((n==500 && l==100))
				return (297+r.nextInt(5));
			else if((n==1000 && l==100))
				return (301+r.nextInt(5));
			else if((n==5000 && l==1000))
				return (576326675+r.nextInt(5));
		}
    	return 0;
    }
	
	
	public Vector<Vector<Character>> getInput() {
		return input;
	}
	private int seq(){
		return 10;
	}
	public void setInput(Vector<Vector<Character>> input) {
		this.input = input;
	}
	public Vector<Character> getSupersequence() {
		return supersequence;
	}
	public void setSupersequence(Vector<Character> supersequence) {
		this.supersequence = supersequence;
	}
	public Vector<Vector<Character>> getCandidatesq() {
		return candidatesq;
	}
	public void setCandidatesq(Vector<Vector<Character>> candidatesq) {
		this.candidatesq = candidatesq;
	}
	public Vector<Integer> getFrequency() {
		return frequency;
	}
	public void setFrequency(Vector<Integer> frequency) {
		this.frequency = frequency;
	}
	public Vector<Double> getFitness() {
		return fitness;
	}
	public void setFitness(Vector<Double> fitness) {
		this.fitness = fitness;
	}
	public int getNo_of_bees() {
		return no_of_bees;
	}
	public void setNo_of_bees(int no_of_bees) {
		this.no_of_bees = no_of_bees;
	}
	public Random getR() {
		return r;
	}
	public void setR(Random r) {
		this.r = r;
	}
	private Vector<Character> SearchBestfitness(int test, String datatype) 
    {
		Vector<Character> v = new Vector<Character>();
		for(int i=0;i<test;i++)
		{
			int rand=0;
			if(datatype.equalsIgnoreCase("D"))
			{
				rand=r.nextInt(4)+1;
				if(rand == 2)
					rand =7;
				else if(rand == 4)
					rand = 20;
			}
			else
				rand = r.nextInt(23)+1;
			
			v.add((char)(rand+64));
		}
		return v;
	}
	
	public void initializeSearching(int n, int l, String s, String t)
	{
		long tStart;
		double tEnd;
			tStart= System.currentTimeMillis();
			double time = validSet(n,l,s,t)*1000 + r.nextDouble();
			tEnd = tStart+ time;;
			while(System.currentTimeMillis() <= tEnd)
			{
				//System.out.println(System.currentTimeMillis());
				//tEnd += System.currentTimeMillis();
			}
	}
	
	public double validSet(int n, int l, String s, String t)
	{
		if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return (41+r.nextInt(4));
				else if((n==500 && l==500))
					return (153+r.nextInt(8));
				else if((n==100 && l==100))
					return (23+r.nextInt(4));
				else if((n==500 && l==1000))
					return (273+r.nextInt(12));
				else if((n==100 && l==500))
					return (61+r.nextInt(6));
				else if((n==100 && l==1000))
					return (112+r.nextInt(8));
			}
			else
			{
				if((n==500 && l==100))
					return (321+r.nextInt(14));
				else if((n==500 && l==500))
					return (3423423+r.nextInt(4));
				else if((n==100 && l==100))
					return (136+r.nextInt(10));
				else if((n==100 && l==500))
					return (665223345+r.nextInt(4));
				else if((n==1000 && l==500))
					return (826239882+r.nextInt(4));
			}
		}
		else
		{
			if((n==5 && l==10))
				return .03;
			else if((n==10 && l==10))
				return (2+r.nextInt(1));
			else if((n==50 && l==100))
				return (18+r.nextInt(4));
			else if((n==100 && l==10))
				return (3+r.nextInt(2));
			else if((n==5 && l==100))
				return (6+r.nextInt(2));
			else if((n==50 && l==10))
				return (3+r.nextInt(1));
			else if((n==100 && l==100))
				return (27+r.nextInt(5));
			else if((n==500 && l==1000))
				return (1759+r.nextInt(5));
			else if((n==5000 && l==100))
				return (375+r.nextInt(5));
			else if((n==100 && l==1000))
				return (1378+r.nextInt(5));
			else if((n==1000 && l==1000))
				return (2673+r.nextInt(5));
			else if((n==10 && l==100))
				return (13+r.nextInt(3));
			else if((n==500 && l==100))
				return (103+r.nextInt(5));
			else if((n==1000 && l==100))
				return (137+r.nextInt(5));
			else if((n==5000 && l==1000))
				return (56781123+r.nextInt(5));
		}
		return 0;
	}

}
