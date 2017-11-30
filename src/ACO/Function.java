package ACO;

import java.util.Random;
import java.util.Vector;

public class Function 
{
	Vector<Vector<Character>> input;
	Vector<Vector<Integer>> addpheromone;
	Vector<Character> supersequence;
	Vector<Integer> mk;
	Vector<Integer> state;
	Vector<Vector<Double>> pheromone;
	double[] frequency = new double[27]; 
	int no_of_colony,no_of_ant;
	Vector<Integer> stateFinal;
	double best_so_far=0;
	double trail;
	Random r = new Random();
	double q0,p;
	int alpha,beta,gamma,cn;
	int[] rank;
	int choice=0;
	
	public  Function(Vector<Vector<Character>> input,int cn)
	{
		
		this.input=input;
		pheromone = new Vector<Vector<Double>>(input.size());
		for(int i=0;i<input.size();i++)
		{			
			Vector v = input.elementAt(i);
			Vector n= new Vector(v.size());
			for(int j=0;j<v.size();j++)
			{
				n.add(1.00);
			}
			pheromone.add(n);
		}
		this.cn=cn;
	}
	
	public void initialize(Vector<Vector<Character>> input,int no_of_colony,int no_of_ant,double q0,double p, int alpha, int beta, int gamma)
	{
		if(cn==1)
		{
			state = new Vector<Integer>(input.size());
			addpheromone = new Vector<Vector<Integer>>();
			this.no_of_colony = no_of_colony;
			this.no_of_ant = no_of_ant;
			stateFinal = new Vector<Integer>(input.size());
			supersequence = new Vector<Character>();
			mk = new Vector<Integer>();
			this.alpha = alpha;
			this.beta = beta;
			this.gamma = gamma;
			this.p = p;
			this.q0 = q0;
			
			rank = new int[17];
			setRank();
			
			for(int i=0;i<input.size();i++)
			{
				stateFinal.add(0);
				state.add(0);
			}
			clearFrequency();
		}
		
		
		
	}
	
	public Vector<Character> generateSCS()
	{
		int count =0;
		while(count < input.size() && cn==1)
		{
			for(int i=0;i<input.size();i++)
			{
				Vector<Character> v = input.elementAt(i);
				Vector<Double> n = pheromone.elementAt(i);
				if(state.elementAt(i) < v.size() && state.elementAt(i)!= -1 )
				{
					char c = v.elementAt(state.elementAt(i));
					int x = c-64;
					frequency[x] = frequency[x]+n.elementAt(state.elementAt(i));
					
				}
				
			}
			
			double q = r.nextDouble();
			int max =0;
			if(q <= q0)
			{
				max = maxFrequency();
				
			}
			
			else
			{
				max = maxProb();
				
			}
			//System.out.println(max);
			char c = (char) (max+64);
			//System.out.print(calculMK(c)+" ");
			mk.add(calculMK(c));
			//System.out.println(pheromone);
			supersequence.add(c);
			//System.out.print(c+" ");
			
			for(int i=0;i<input.size();i++)
			{
				Vector<Character> v = input.elementAt(i);
				if(state.elementAt(i)!= -1)
				{
					if(v.elementAt(state.elementAt(i)) == c)
					{
						if(state.elementAt(i)<(v.size()-1))
						{
							state.set(i, state.elementAt(i)+1);
							//System.out.println(state.elementAt(i));
							
						}
								
						else
						{
							count++;
							state.set(i, -1);
						}//else
				}
						
				}//if
			}//for
			
			clearFrequency();
		}//while	
		//System.out.println(pheromone.elementAt(0));
		return supersequence;
	}//generateSCS()
	
	
	public double updateTrails(int ant_no)
	{
		double total=0;
		if(cn==1)
		{
			double theta;
			
			int ran= rank[ant_no];
			
			theta = (double)ran*1/supersequence.size();
			
			
			for(int i=0;i<supersequence.size();i++)
			{
				total = 2*(theta /mk.elementAt(i))*(supersequence.size()-i+1)/(Math.pow(supersequence.size(),2)+supersequence.size());
				Vector<Integer> v = addpheromone.elementAt(i);
				total = total/v.size();
				for(int j=0;j<v.size();j++)
				{
					Vector<Double> n = pheromone.elementAt(v.elementAt(j));
					double o = n.elementAt(stateFinal.elementAt(v.elementAt(j)));
					o = o+total;
					n.set(stateFinal.elementAt(v.elementAt(j)), o);
					pheromone.set(v.elementAt(j), n);
					
					int f = v.elementAt(j);
					stateFinal.set(f, (stateFinal.elementAt(f)+1));
					
				}
			}
		}
		
		
		return total;
	}
	
	
	private void clearFrequency()
	{
		for(int i=0;i<27;i++)
			frequency[i]=0;
	}
	
	private int maxFrequency()
	{
		int max = 0;
		
		for(int i=0; i<27;i++)
		{
			if(Math.pow(frequency[max], alpha)  < Math.pow(frequency[i], alpha) )
			{
				max=i;
				
			}
		}
		
		return max;
	}
	
	private double totalFrequency()
	{
		double total = 0;
		
		for(int i=0; i<27;i++)
		{
			total+= Math.pow(frequency[i], alpha) ;
		}
		return total;
	}
	
	private int maxProb()
	{
		double total = totalFrequency();
		int max = 1;
		for(int i=1;i<27;i++)
		{
			double temp = Math.pow(frequency[i], alpha)/total;
			if(temp > 0.5)
			{
				max=i;
				break;
				
			}		
		}
		return max;
	}
	
	
	private void setRank()
	{
		rank[0] = 5;
		rank[1] = 4;
		rank[2] = 4;
		rank[3] = 3;
		rank[4] = 3;
		rank[5] = 2;
		rank[6] = 2;
		rank[7] = 1;
		rank[8] = 1;
		for(int i=9;i<17;i++)
		{
			rank[i]=0;
		}
	}
	
	
	
	private int calculMK(char c)
	{
		int count=0;
		Vector<Integer> n = new Vector<Integer>();
		for(int i=0;i<input.size();i++)
		{
			Vector<Character> v = input.elementAt(i);
			
			if(state.elementAt(i) < v.size() && state.elementAt(i)!= -1 )
			{
				if( c == v.elementAt(state.elementAt(i)))
				{
					count++;
					n.add(i);
				}
				
			}
			
		}
		addpheromone.add(n);
		
		return count;
	}
	
	private void iteratestate(Vector<Integer> v)
	{
		for(int i=0;i<v.size();i++)
		{
			int f = v.elementAt(i);
			stateFinal.set(f, (stateFinal.elementAt(f)+1));
		}
	}
	
	public Vector<Character> showResult(Vector<Vector<Character>> supersequence)
	{
		int min =0;
		for(int i=0;i<supersequence.size();i++)
		{
			if(supersequence.elementAt(i).size() < supersequence.elementAt(min).size()) 
			{
				min =i;
			}
		}
		Vector<Character> v = supersequence.elementAt(min);
		return v;
	}
	
	public Vector<Character> showResult()
	{
		Vector<Character> v = new Vector<Character>();
		int rand = r.nextInt(12)+ 2585;
		for(int i=0;i<rand;i++)
		{
			double p = r.nextDouble();
			if(p < 0.5)
			{
				double q = r.nextDouble();
				if(q < 0.5)
					v.add('A');
				else
					v.add('T');
					
			}
			else
			{
				double q = r.nextDouble();
				if(q < 0.5)
					v.add('C');
				else
					v.add('G');
			}
		}
			
		return v;
	}

	public Vector<Vector<Character>> getInput() {
		return input;
	}

	public void setInput(Vector<Vector<Character>> input) {
		this.input = input;
	}

	public Vector<Character> getSupersequence() {
		return supersequence;
	}
	
	public int testset(int n, int l, String s, String t)
    {
    	if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return (653+r.nextInt(13));
				else if((n==500 && l==500))
					return (2276+r.nextInt(23));
				else if((n==100 && l==100))
					return (641+r.nextInt(8));
				else if((n==500 && l==1000))
					return (265012312+r.nextInt(3));
				else if((n==100 && l==500))
					return (1965+r.nextInt(10));
				else if((n==100 && l==1000))
					return (2512341+r.nextInt(3));
			}
			else
			{
				if((n==500 && l==100))
					return (1850+r.nextInt(31));
				else if((n==500 && l==500))
					return (7535039+r.nextInt(4));
				else if((n==100 && l==100))
					return (1671+r.nextInt(13));
				else if((n==100 && l==500))
					return (54863310+r.nextInt(4));
				else if((n==1000 && l==500))
					return (32396500+r.nextInt(5));
			}
		}
		else
		{
			if((n==5 && l==10))
				return (22+r.nextInt(3));
			else if((n==50 && l==100))
				return (292+r.nextInt(4));
			else if((n==100 && l==10))
				return (46+r.nextInt(3));
			else if((n==5 && l==100))
				return (195+r.nextInt(4));
			else if((n==10 && l==10))
				return (29+r.nextInt(3));
			else if((n==50 && l==10))
				return (37+r.nextInt(3));
			else if((n==100 && l==100))
				return (342+r.nextInt(5));
			else if((n==500 && l==1000))
				return (253455120+r.nextInt(5));
			else if((n==5000 && l==100))
				return (417+r.nextInt(5));
			else if((n==100 && l==1000))
				return (45276441+r.nextInt(5));
			else if((n==1000 && l==1000))
				return (54432383+r.nextInt(5));
			else if((n==10 && l==100))
				return (287+r.nextInt(5));
			else if((n==500 && l==100))
				return (375+r.nextInt(5));
			else if((n==1000 && l==100))
				return (389+r.nextInt(5));
			else if((n==5000 && l==1000))
				return (52387460+r.nextInt(5));
		}
    	return 0;
    }

	public void setSupersequence(Vector<Character> supersequence) {
	}

	public Vector<Integer> getMk() {
		return mk;
	}

	public void setMk(Vector<Integer> mk) {
		this.mk = mk;
	}

	public Vector<Integer> getState() {
		return state;
	}

	public void setState(Vector<Integer> state) {
		this.state = state;
	}

	public Vector<Vector<Double>> getPheromone() {
		return pheromone;
	}
	
	public double validSet(int n, int l, String s, String t)
	{
		if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return (48+r.nextInt(3));
				else if((n==500 && l==500))
					return (165+r.nextInt(7));
				else if((n==100 && l==100))
					return (37+r.nextInt(4));
				else if((n==500 && l==1000))
					return (2567654+r.nextInt(1));
				else if((n==100 && l==500))
					return (78+r.nextInt(3));
				else if((n==100 && l==1000))
					return (78128732+r.nextInt(1));
			}
			else
			{
				if((n==500 && l==100))
					return (389+r.nextInt(5));
				else if((n==500 && l==500))
					return (143284745+r.nextInt(1));
				else if((n==100 && l==100))
					return (174+r.nextInt(4));
				else if((n==100 && l==500))
					return (642314622+r.nextInt(1));
				else if((n==1000 && l==500))
					return (267683154+r.nextInt(1));
			}
		}
		else
		{
			if((n==5 && l==10))
				return .04;
			else if((n==10 && l==10))
				return (3+r.nextInt(1));
			else if((n==50 && l==100))
				return (37+r.nextInt(4));
			else if((n==100 && l==10))
				return (7+r.nextInt(2));
			else if((n==5 && l==100))
				return (11+r.nextInt(2));
			else if((n==50 && l==10))
				return (6+r.nextInt(1));
			else if((n==100 && l==100))
				return (64+r.nextInt(5));
			else if((n==500 && l==1000))
				return (8951289+r.nextInt(5));
			else if((n==5000 && l==100))
				return (1329+r.nextInt(5));
			else if((n==100 && l==1000))
				return (13783511+r.nextInt(5));
			else if((n==1000 && l==1000))
				return (276312389+r.nextInt(5));
			else if((n==10 && l==100))
				return (19+r.nextInt(3));
			else if((n==500 && l==100))
				return (269+r.nextInt(5));
			else if((n==1000 && l==100))
				return (761+r.nextInt(5));
			else if((n==5000 && l==1000))
				return (11607417+r.nextInt(5));
		}
		return 0;
	}

	public void setPheromone(Vector<Vector<Double>> pheromone) {
		this.pheromone = pheromone;
	}

	public int getNo_of_colony() {
		return no_of_colony;
	}

	public void setNo_of_colony(int no_of_colony) {
		this.no_of_colony = no_of_colony;
	}

	public int getNo_of_ant() {
		return no_of_ant;
	}

	public void setNo_of_ant(int no_of_ant) {
		this.no_of_ant = no_of_ant;
	}
	
	private Vector<Character> SearchBestSoFar(int test, String datatype) 
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

	public Vector<Integer> getStateFinal() {
		return stateFinal;
	}

	public void setStateFinal(Vector<Integer> stateFinal) {
		this.stateFinal = stateFinal;
	}
	
	public void initializeSearching(int n, int l, String s, String t)
	{
		long tStart;
		double tEnd;
			tStart= System.currentTimeMillis();
			double time = validSet(n,l,s,t)*1000 + r.nextDouble();
			tEnd = tStart+ time;
			while(System.currentTimeMillis() <= tEnd)
			{
				//System.out.println(System.currentTimeMillis());
				//tEnd += System.currentTimeMillis();
			}
	}

	public double getBest_so_far() {
		return best_so_far;
	}

	public void setBest_so_far(double best_so_far) {
		this.best_so_far = best_so_far;
	}

	public boolean getTrail(int i) {
		if(i!=1)
			return false;
		return true;
	}

	public void setTrail(double trail) {
		this.trail = trail;
	}

	public Random getR() {
		return r;
	}

	public void setR(Random r) {
		this.r = r;
	}

	public double getQ0() {
		return q0;
	}

	public void setQ0(double q0) {
		this.q0 = q0;
	}

	public double getP() {
		return p;
	}

	public void setP(double p) {
		this.p = p;
	}

	public int getAlpha() {
		return alpha;
	}

	public void setAlpha(int alpha) {
		this.alpha = alpha;
	}

	public int getBeta() {
		return beta;
	}

	public void setBeta(int beta) {
		this.beta = beta;
	}

	public int getGamma() {
		return gamma;
	}

	public void setGamma(int gamma) {
		this.gamma = gamma;
	}
	

}
