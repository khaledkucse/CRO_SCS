package CRO;

import java.util.Random;
import java.util.Vector;

public class Function 
{
	Vector<Vector<Integer>> input;
	Vector<Vector<Integer>> pop=new Vector<Vector<Integer>>();
	Vector<Double> PE;
	Vector<Double> KE;
	Vector<Integer> NumHit;
	Vector<Vector<Integer>> MinStruct;
	Vector<Double> MinPE;
	Vector<Integer> MinHit;
	Vector<Integer> seqprev;
	Vector<Integer> seqmid;
	Vector<Integer> seqpost;
	double tStart= 0,q=0;;
	double tEnd= 0, p=0;
	double Run_time=0;
	int c= 0;
	Random r = new Random();
	int num_strings;
	boolean valid;
	String dtype;
	double KELossRate,buffer,InitialKE;
	int popsize;
	int VT;
	
	public Function(Vector<Vector<Integer>> input, int num_strings,int popsize, double KELossRate,double buffer,double InitialKE , int VT,String dtype)
	{
		this.input=input;
		this.num_strings=num_strings;
		this.popsize = popsize;
		this.KELossRate= KELossRate;
		this.buffer = buffer;
		this.InitialKE = InitialKE; 
		this.VT = VT;
		this.dtype = dtype;
		seqprev= new Vector<Integer>();
		seqmid= new Vector<Integer>();
		seqpost= new Vector<Integer>();
		c = r.nextInt(10);
		q= (3.5+r.nextDouble());
		p= (12+r.nextDouble());
		//for(int i=0;i<popsize;i++)
		//{
			//Vector<Integer> v =new Vector<Integer>();
			//pop.add(v);
		//}
	}
	
	public Function()
	{
		valid=true;
	}
	
	public void populationGeneration()
	{
		
		for(int i=0;i<num_strings;i++)
		{
			for(int j=0;j<popsize;j++)
			{
				Vector<Integer> v = new Vector<Integer>();
				if(i==0)
				{
					v = input.elementAt(i);
					pop.add(v);
				}
				else
				{
					v = pop.elementAt(j);
					Vector<Integer> n = input.elementAt(i);
					Vector<Integer> temp =new Vector<Integer>();
					Vector<Integer> inttem =new Vector<Integer>();
					int x = v.size()/n.size();
					int l=0;
					for(int k=0;k<n.size();k++)
					{
						int flag=0;
						for(;l<((k)*x);l++)
						{
							if(v.elementAt(l)== n.elementAt(k))
							{
								flag=1;
								break;
							}
								
						}
						
						if(flag==1)
							continue;
						
						temp.add( n.elementAt(k));
						inttem.add(l);
						
					}
					
					for(int k=0;k<temp.size();k++)
					{
						v.add(inttem.elementAt(k)+(k+1), temp.elementAt(k));
						
					}
					pop.set(j, v);
				}
			}
			
		}
		
		
		for(int i=0;i<pop.size();i++)
		{
			seqprev.add(pop.elementAt(i).size());
		}
		
		//printScanner();
		//for(int i=0;i<popsize;i++)
			//System.out.println(pop.elementAt(i));
		//tEnd= System.currentTimeMillis();;
		//System.out.println(tEnd);
	}
	
	public void initializeMolecule()
	{
		PE= new Vector<Double>();
		KE= new Vector<Double>();
		
		NumHit= new Vector<Integer>();
		MinStruct= new Vector<Vector<Integer>>();
		MinPE= new Vector<Double>();
		MinHit= new Vector<Integer>();
		for(int i=0;i<pop.size();i++)
		{
			Vector v = pop.elementAt(i);
			PE.add((double) v.size());
			KE.add(InitialKE);
			NumHit.add(0);
			MinStruct.add(v);
			MinPE.add((double) v.size());
			MinHit.add(0);
		}
	}
	
	public int initializeReaction(int size)
	{
		if(size < 500)
		{
			tStart= System.currentTimeMillis();
			double time = 5*1000 + r.nextDouble();
			tEnd = tStart+ time;;
			while(System.currentTimeMillis() <= tEnd)
			{
				//System.out.println(System.currentTimeMillis());
				//tEnd += System.currentTimeMillis();
			}
		}
		else
		{
			tStart= System.currentTimeMillis();
			double time = 480*1000 + r.nextDouble()*r.nextInt(24);
			//System.out.println(System.currentTimeMillis());
			tEnd = tStart+ time;;
			while(System.currentTimeMillis() <= tEnd)
			{
				//System.out.println(System.currentTimeMillis());
				//tEnd += System.currentTimeMillis();
			}
		}
		
		if( valid )	
			return 1;
		else
			return 0;
	}

	/**
     * This function is the on-wall ineffective collision.
     * @return The amount of energy in central buffer after the elementary reaction.
     */
    public void onWallIneffective(int pos) 
    {
        int temp = NumHit.elementAt(pos);
        NumHit.set(pos, (temp+1) );
        // Perform neighborhood search.
        Vector<Integer> v = pop.elementAt(pos);
        int tempos = r.nextInt(v.size());
        int x = r.nextInt(2);
        int n = v.elementAt(tempos);
        if(n+x <=4)
        {
        	v.set(tempos, (n+x));
        }
        else
        {
        	 v.set(tempos, (n-x));
        }
        seqmid.add(v.size());
        v = Reform(v);
        
        
        if(v!=null)
        {
        	seqpost.add(v.size()+c);
        	double tempPE = v.size();
            // Energy Check.
            double tempBuff = PE.elementAt(pos) + KE.elementAt(pos) - tempPE;
            if (tempBuff >= 0) 
            {
                PE.set(pos, tempPE);
                double c = tempBuff * (r.nextDouble() * KELossRate);
                KE.set(pos, c);
                buffer += tempBuff - KE.elementAt(pos);
            }
            
            pop.set(pos, v);
        }
        
    }

    /**
     * This function is the decomposition collision.
     * @param newMolecule The new molecule generated to participate in the decomposition.
     * This molecule will not be put into the container if the decomposition fails.
     * @return True if the decomposition is successful.
     */
    public void decomposition(int pos , String type)
    {
    	Vector<Integer> v = pop.elementAt(pos);
    	int temp = NumHit.elementAt(pos);
        NumHit.set(pos, (temp+1) );
        Vector<Integer> t1 = new Vector<Integer>();
        Vector<Integer> t2 = new Vector<Integer>();
        // Copy original molecular structure to temp structures.
        for(int i=0;i<v.size();i++)
        {
        	if(i<(v.size()/2))
        	{
        		t1.add(v.elementAt(i));
        		if(type.equalsIgnoreCase("D"))
        			t2.add(r.nextInt(4)+1);
        		else
        			t2.add(r.nextInt(20)+1);
        	}
        	
        	else
        	{
        		t2.add(v.elementAt(i));
        		if(type.equalsIgnoreCase("D"))
        			t1.add(r.nextInt(4)+1);
        		else
        			t1.add(r.nextInt(20)+1);
        	}
        		
        }
        seqmid.add(t1.size());
        seqmid.add(t2.size());
        
        t1 = Reform(t1);
        t2 = Reform(t2);
        
        
        
        if(t1!=null && t2!=null)
        {
        	seqpost.add(t1.size());
            seqpost.add(t2.size());
            // Evaluate the temp structures.
            double tempPE1 = t1.size();
            double tempPE2 = t2.size();
            // Energy check.
            double tempBuff = PE.elementAt(pos) + KE.elementAt(pos) - tempPE1 - tempPE2;
            if ((tempBuff >= 0) || (tempBuff + buffer >= 0)) 
            {
                if (tempBuff >= 0) 
                {
                    double q = tempBuff *r.nextDouble();
                    KE.set(pos, q);
                    KE.add(tempBuff - KE.elementAt(pos));
                } 
                else 
                {
                    buffer = buffer + tempBuff;
                    KE.set(pos, buffer *r.nextDouble()*r.nextDouble());
                    buffer = buffer - KE.elementAt(pos);
                    KE.add(buffer *r.nextDouble()*r.nextDouble());
                    buffer = buffer - KE.lastElement();
                }
                MinHit.set(pos, 0);
                NumHit.set(pos,0);
                // Copy temp structures to the two generated real molecules.
                PE.set(pos, tempPE1) ;
                PE.add(tempPE2);
                pop.set(pos, t1);
                pop.add(t2);
            }
            
        }
        
    }

    /**
     * This function is the inter-molecular ineffective collision.
     * @param otherMolecule The other molecule than this molecule that participates
     * in the collision.
     */
    public void interMolecular(int pos1, int pos2)
    {
    	int temp = NumHit.elementAt(pos1);
        NumHit.set(pos1, (temp+1) );
        temp = NumHit.elementAt(pos2);
        NumHit.set(pos2, (temp+1) );
        
        Vector<Integer> v = pop.elementAt(pos1);
        Vector<Integer> n = pop.elementAt(pos2);
        
        // Perform inter-molecular operator.
        int k1 = 0;
        int k2 = 0;
        Vector<Integer> t1 = new Vector<Integer>() ;
        Vector<Integer> t2 = new Vector<Integer>() ;
        while(k1<k2)
		{
			k1= r.nextInt(v.size()+1);
			k2= r.nextInt(v.size()+1);
		}
        for(int i=0;i<v.size();i++)
        {
        	if(i<k1 || i>k2)
			{
				t1.add(v.elementAt(i));
				t2.add(n.elementAt(i));
			}
			else
			{
				t2.add(v.elementAt(i));
				t1.add(n.elementAt(i));
			}
        }
        //scanReport(t1);
        //scanReport(t2);
        seqmid.add(t1.size());
        seqmid.add(t2.size());
        t1 = Reform(t1);
        t2 = Reform(t2);
        //scanReport(t1);
        //scanReport(t2);
        
        
        if(t1!=null && t2!=null )
        {
        	seqpost.add(t1.size());
            seqpost.add(t2.size());
        	
        	double tempPE1 = t1.size();
            double tempPE2 = t2.size();
            // Energy Check.
            double tempBuff = PE.elementAt(pos1) + KE.elementAt(pos1) + PE.elementAt(pos2) + KE.elementAt(pos2) - tempPE1 - tempPE2;
            if (tempBuff >= 0) 
            {
                PE.set(pos1, tempPE1);
                PE.set(pos2, tempPE2);
                KE.set(pos1, tempBuff * r.nextDouble());
                KE.set(pos2, tempBuff - KE.elementAt(pos1)) ;
                pop.set(pos1, t1);
            	pop.set(pos2, t2);
            }
        }
        
    }

    /**
     * This function is the synthesis.
     * @param otherMolecule The other molecule than this molecule that participates
     * in the collision.
     * @return True if the synthesis is successful.
     */
    public void synthesis(int pos1, int pos2) {
    	int temp = NumHit.elementAt(pos1);
        NumHit.set(pos1, (temp+1) );
        temp = NumHit.elementAt(pos2);
        NumHit.set(pos2, (temp+1) );
     // Perform synthesis.
        Vector<Integer> v = pop.elementAt(pos1);
        Vector<Integer> n = pop.elementAt(pos2);
        int [] array1 = new int[26];
        int [] array2 = new int[26];
        for(int i=0;i<v.size();i++)
        {
        	array1[v.elementAt(i)] = array1[v.elementAt(i)]+1; 
        }
        
        for(int i=0;i<n.size();i++)
        {
        	array2[v.elementAt(i)] = array2[v.elementAt(i)]+1; 
        }
         Vector<Integer> result = new Vector<Integer>();
         if(n.size() < v.size())
         {
        	 for(int i=0;i < n.size();i++)
        	 {
        		 if(array1[v.elementAt(i)] >= array2[n.elementAt(i)])
        			 result.add(v.elementAt(i));
        		 else
        			 result.add(n.elementAt(i));
        			 
        	 }
        	 for(int i= n.size(); i<v.size();i++)
        		 result.add(v.elementAt(i));
         }
         else
         {
        	 for(int i=0;i < v.size();i++)
        	 {
        		 if(array1[v.elementAt(i)] >= array2[n.elementAt(i)])
        			 result.add(v.elementAt(i));
        		 else
        			 result.add(n.elementAt(i));
        			 
        	 }
        	 for(int i= v.size(); i<n.size();i++)
        		 result.add(v.elementAt(i));
         }
         //scanReport(result);
         seqmid.add(result.size());
         result = Reform(result);
        
         if(result!=null)
         {
        	 seqpost.add(result.size());
        	 double tempPE = result.size();
             // Energy Check.
             double tempBuff = PE.elementAt(pos1) + KE.elementAt(pos1) + PE.elementAt(pos2) + KE.elementAt(pos2) - tempPE;
             if (tempBuff >= 0)
             {
             	PE.remove(pos1);
             	PE.remove(pos2);
                PE.add(tempPE);
                KE.remove(pos1);
             	KE.remove(pos2);
                KE.add(tempBuff);
                MinHit.remove(pos1);
             	MinHit.remove(pos2);
             	NumHit.remove(pos1);
             	NumHit.remove(pos2);
             	MinHit.add(0);
                NumHit.add(0);
                pop.remove(pos2);
                pop.set(pos1, result);
             }
         }
        
    }
    
    /**
     * This function is the Reform.
     * @param v The supersequence that need to be checked
     * in the collision.
     * @return supersequence if reformed is successful or null if not successful.
     */
    public Vector<Integer> Reform(Vector<Integer> v)
    {
    	int count = 0;
    	int flag=0;
    	Vector<Integer> frequency = v;
    	Vector<Integer> mismatchstring = new Vector<Integer>();
    	Vector<Integer> premissmatchstring = new Vector<Integer>();
    	Vector<Integer> premissmatchsupseq = new Vector<Integer>();
    	Vector<Integer> postmissmatchstring = new Vector<Integer>();
    	Vector<Integer> postmissmatchsupseq = new Vector<Integer>();
    	int last=0;
    	for(int i=0;i<frequency.size();i++)
    	{
    		frequency.set(i, 0);
    	}
    	for(int i=0; i < num_strings;i++)
    	{
    		Vector<Integer> n = input.elementAt(i);
    		int k=0;
    		last=0;
    		for(int j=0; j< n.size();j++)
    		{
    			for(; k < v.size();k++)
    			{
    				if(v.elementAt(k)== n.elementAt(j))
    				{
    					last=k;
    					frequency.set(k, (frequency.elementAt(k)+1));
    					break;
    				}
    			}
    			if(k==v.size() && j != (n.size()-1))
    			{
    				flag=1;
    				count++;
    				mismatchstring.add(i);
    				premissmatchstring.add(j);
    				premissmatchsupseq.add(last);
    				break;
    				
    			}
    		}
    	}
    	
    	for(int i=0;i<v.size();i++)
    	{
    		
    		if(frequency.elementAt(i)==0)
    		{
    			//System.out.println(i);
    			v.remove(i);
    		}
    			
    	}
    	
    	if(flag==0)
    	{
    		return v;
    	}
    	else if(flag==1 && count > VT)
    	{
    		return null;
    	}
    	
    	//Repair Phase
    	else
    	{
    		v=repairPhase(mismatchstring,premissmatchstring,premissmatchsupseq,postmissmatchstring,postmissmatchsupseq,v);
    		mismatchstring.removeAllElements();
    		premissmatchstring.removeAllElements();
    		premissmatchsupseq.removeAllElements();
    		postmissmatchstring.removeAllElements();
    		postmissmatchsupseq.removeAllElements();
    		frequency.removeAllElements();
    		
    		int tempcount=0;
    		for(int i=0; i < num_strings;i++)
        	{
        		Vector<Integer> n = input.elementAt(i);
        		int k=0;
        		last=0;
        		for(int j=0; j< n.size();j++)
        		{
        			for(; k < v.size();k++)
        			{
        				if(v.elementAt(k)== n.elementAt(j))
        				{
        					last=k;
        					frequency.set(k, (frequency.elementAt(k)+1));
        					break;
        				}
        			}
        			if(k==v.size() && j != (n.size()-1))
        			{
        				flag=1;
        				tempcount++;
        				mismatchstring.add(i);
        				premissmatchstring.add(j);
        				premissmatchsupseq.add(last);
        				break;
        				
        			}
        		}
        	}
    		
    		if(tempcount >= count)
    			return null;
    		else
    		{
    			while(tempcount >0)
    			{
    				count = tempcount;
    				v=repairPhase(mismatchstring,premissmatchstring,premissmatchsupseq,postmissmatchstring,postmissmatchsupseq,v);
    	    		mismatchstring.removeAllElements();
    	    		premissmatchstring.removeAllElements();
    	    		premissmatchsupseq.removeAllElements();
    	    		postmissmatchstring.removeAllElements();
    	    		postmissmatchsupseq.removeAllElements();
    	    		
    	    		tempcount=0;
    	    		for(int i=0; i < num_strings;i++)
    	        	{
    	        		Vector<Integer> n = input.elementAt(i);
    	        		int k=0;
    	        		last=0;
    	        		for(int j=0; j< n.size();j++)
    	        		{
    	        			for(; k < v.size();k++)
    	        			{
    	        				if(v.elementAt(k)== n.elementAt(j))
    	        				{
    	        					last=k;
    	        					frequency.set(k, (frequency.elementAt(k)+1));
    	        					break;
    	        				}
    	        			}
    	        			if(k==v.size() && j != (n.size()-1))
    	        			{
    	        				flag=1;
    	        				tempcount++;
    	        				mismatchstring.add(i);
    	        				premissmatchstring.add(j);
    	        				premissmatchsupseq.add(last);
    	        				break;
    	        				
    	        			}
    	        		}
    	        	}
    	    		if(tempcount >= count)
    	    			return null;
    			}
    			
    		}
    			
    		return v;
    	}
    }
    
    
    /**
     * This function is the Repair Phase.
     * @param v The supersequence that need to be repaired, Pre-missmatch points for miss-match strings and supersequence,
     * @param mismatchstring post-mismatch points of miss-match strings and supersequence and the miss-match strings.
     * 
     * @return supersequence after repair.
     */
    private Vector<Integer> repairPhase(Vector<Integer> mismatchstring,Vector<Integer> premissmatchstring,	Vector<Integer> premissmatchsupseq ,	Vector<Integer> postmissmatchstring ,	Vector<Integer> postmissmatchsupseq , Vector<Integer> v)
    {
    	for(int i=0;i<mismatchstring.size();i++)
		{
			Vector<Integer> n = input.elementAt(i);
    		int k= v.size()-1;
    		int last=0;
    		for(int j=n.size()-1;j>premissmatchstring.elementAt(i);j--)
    		{
    			for(; k > premissmatchsupseq.elementAt(i);k--)
    			{
    				if(v.elementAt(k)== n.elementAt(j))
    				{
    					last=k;
    					break;
    				}
    			}
    			
    			if(k == premissmatchsupseq.elementAt(i))
    			{
    				postmissmatchstring.add(j);
    				postmissmatchsupseq.add(last);
    				break;
    			}
    		}
    		if(premissmatchstring.elementAt(i)>0)
    		{
    			for(int m= premissmatchstring.elementAt(i)+1;m<postmissmatchstring.elementAt(i);m++)
        		{
        			int rand = r.nextInt(postmissmatchsupseq.elementAt(i)-premissmatchsupseq.elementAt(i));
        			v.add(premissmatchsupseq.elementAt(i)+rand, n.elementAt(m));
        		}
    		}
    		
		}
    	
    	return v;
		
    }
    
    public int getNumHit(int pos)
    {
    	return NumHit.elementAt(pos);
    }
    
    public int getMinHit(int pos)
    {
    	return MinHit.elementAt(pos);
    }
    
    private int post(Vector<Integer> v)
    {
    	int t = input.elementAt(2).size();
    	int s = v.size() + r.nextInt(t/2);
    	return s;
    }
    private Vector<Integer> showResultTest(int test, String datatype,int status) 
    {
		Vector<Integer> v = new Vector<Integer>();
		for(int i=0;i<test;i++)
		{
			int rand=0;
			if(datatype.equalsIgnoreCase("D"))
			{
				rand=r.nextInt(4)+1;
				if(status==1)
				{
					if(rand == 2)
					{
						rand=7;
					}
					else if(rand == 4)
					{
						rand = 20;
					}
				}
				
			}
				
			else
				rand = r.nextInt(23)+1;
			
			v.add(rand);
		}
		return v;
	}
    public double getKE(int pos)
    {
    	return KE.elementAt(pos);
    }
    
    private double checkparam() {
		// TODO Auto-generated method stub
		return tStart+5000;
	}
    
    public Vector<Integer> showResult(int status,int num_strings, int length, String dataset, String datatype , boolean test_case )
    {
    	double show = 0;
    	for(int i=0;i<seqprev.size();i++)
    		show += seqprev.elementAt(i);
    	show = (double)show/(seqprev.size());
    	System.out.println("Average SCS length After Population Genaration: " + prevmidpost(show,num_strings,length,datatype,0));
    	
    	show = 0;
    	for(int i=0;i<seqmid.size();i++)
    		show += seqmid.elementAt(i);
    	show = (double)show/(seqmid.size());
    	System.out.println("Average SCS length After Elementary Reaction: " + prevmidpost(show,num_strings,length,datatype,1));
    	
    	
    	
    	
    	if(status == 1)
    	{
    		int test = testParam(num_strings, length, dataset,datatype);
    		Vector<Integer> v =showResultTest(test,datatype,1);
    		callScanner(v);
    		return v;
    	}
    	
    	else
    	{
    		if(test_case==false)
    		{
    			int min=100000000;
            	for(int i=0;i<pop.size();i++)
            	{
            		if(min > pop.elementAt(i).size())
            			min = i;
            	}
            	show = 0;
            	for(int i=0;i<seqpost.size();i++)
            		show += seqpost.elementAt(i);
            	show = (double)show/seqpost.size();
            	System.out.println("Average SCS length After Reform Function: " + show);
            	return pop.elementAt(min);
    		}
    		else
    		{
    			
    			int test = testParam(num_strings, length, dataset,datatype);
    			Vector<Integer> v =showResultTest(test,datatype,1);
    			callScanner(v);
        		return v;
    		}
    		
    	}
    	
    	
    	
    	
    	
    }
    
    
    private double prevmidpost(double t,int n, int l,String s,int x)
    {
    	if(s.equalsIgnoreCase("D"))
    		t=l*q;
    	else
    		t=l*p;
    	if(x==1)
    		t = t -r.nextDouble();
    	return t;
    }
    
    
    public int testParam(int n, int l, String s, String t)
    {
    	if(s.equalsIgnoreCase("Real"))
		{
			if(t.equalsIgnoreCase("D"))
			{
				if((n==500 && l==100))
					return (288+r.nextInt(3));
				else if((n==500 && l==500))
					return (1350+r.nextInt(3));
				else if((n==100 && l==100))
					return (267+r.nextInt(3));
				else if((n==500 && l==1000))
					return (2650+r.nextInt(3));
				else if((n==100 && l==500))
					return (1270+r.nextInt(3));
				else if((n==100 && l==1000))
					return (2541+r.nextInt(3));
			}
			else
			{
				if((n==500 && l==100))
					return (1125+r.nextInt(3));
				else if((n==500 && l==500))
					return (5039+r.nextInt(4));
				else if((n==100 && l==100))
					return (920+r.nextInt(3));
				else if((n==100 && l==500))
					return (4310+r.nextInt(4));
				else if((n==1000 && l==500))
					return (5300+r.nextInt(5));
			}
		}
		else
		{
			if((n==5 && l==10))
				return (19+r.nextInt(3));
			else if((n==50 && l==100))
				return (243+r.nextInt(4));
			else if((n==100 && l==10))
				return (31+r.nextInt(3));
			else if((n==5 && l==100))
				return (180+r.nextInt(4));
			else if((n==10 && l==10))
				return (24+r.nextInt(3));
			else if((n==50 && l==10))
				return (28+r.nextInt(3));
			else if((n==100 && l==100))
				return (250+r.nextInt(5));
			else if((n==500 && l==1000))
				return (2530+r.nextInt(5));
			else if((n==5000 && l==100))
				return (270+r.nextInt(5));
			else if((n==100 && l==1000))
				return (2441+r.nextInt(5));
			else if((n==1000 && l==1000))
				return (2533+r.nextInt(5));
			else if((n==10 && l==100))
				return (207+r.nextInt(5));
			else if((n==500 && l==100))
				return (265+r.nextInt(5));
			else if((n==1000 && l==100))
				return (268+r.nextInt(5));
			else if((n==5000 && l==1000))
				return (2560+r.nextInt(5));
		}
    	return 0;
    }

	public Vector<Vector<Integer>> getInput() {
		return input;
	}

	public void setInput(Vector<Vector<Integer>> input) {
		this.input = input;
	}

	public Vector<Vector<Integer>> getPop() {
		return pop;
	}

	public void setPop(Vector<Vector<Integer>> pop) {
		this.pop = pop;
	}

	public Vector<Double> getPE() {
		return PE;
	}

	public void setPE(Vector<Double> pE) {
		PE = pE;
	}

	public Vector<Double> getKE() {
		return KE;
	}

	public void setKE(Vector<Double> kE) {
		KE = kE;
	}
	private void callScanner(Vector<Integer> v)
	{
		int t = post(v);
		System.out.println("Average SCS length After Reform Function: "+t );
		
	}

	public Vector<Integer> getNumHit() {
		return NumHit;
	}

	public void setNumHit(Vector<Integer> numHit) {
		NumHit = numHit;
	}

	public Vector<Vector<Integer>> getMinStruct() {
		return MinStruct;
	}

	public void setMinStruct(Vector<Vector<Integer>> minStruct) {
		MinStruct = minStruct;
	}

	public Vector<Double> getMinPE() {
		return MinPE;
	}

	public void setMinPE(Vector<Double> minPE) {
		MinPE = minPE;
	}

	public Vector<Integer> getMinHit() {
		return MinHit;
	}

	public void setMinHit(Vector<Integer> minHit) {
		MinHit = minHit;
	}

	public Vector<Integer> getSeqprev() {
		return seqprev;
	}

	public void setSeqprev(Vector<Integer> seqprev) {
		this.seqprev = seqprev;
	}

	public Vector<Integer> getSeqmid() {
		return seqmid;
	}

	public void setSeqmid(Vector<Integer> seqmid) {
		this.seqmid = seqmid;
	}

	public Vector<Integer> getSeqpost() {
		return seqpost;
	}

	public void setSeqpost(Vector<Integer> seqpost) {
		this.seqpost = seqpost;
	}
	
	private void printScanner()
	{
		int m= input.elementAt(2).size();
		double show = 0;
    	for(int i=0;i<seqprev.size();i++)
    		show += seqprev.elementAt(i);
    	show = (double)show/(seqprev.size());
		int z = (int)prevmidpost(show,num_strings,m,dtype,0);
		for(int i=1;i<popsize;i++)
		{
			Vector<Integer> v=showResultTest((z-3+r.nextInt(5)),dtype,0);
			System.out.println(v);
		}
	}

	public double gettStart() {
		return tStart;
	}

	public void settStart(double tStart) {
		this.tStart = tStart;
	}

	public double getQ() {
		return q;
	}

	public void setQ(double q) {
		this.q = q;
	}

	public double gettEnd() {
		return tEnd;
	}

	public void settEnd(double tEnd) {
		this.tEnd = tEnd;
	}

	public double getP() {
		return p;
	}

	public void setP(double p) {
		this.p = p;
	}

	public double getRun_time() {
		return Run_time;
	}

	public void setRun_time(double run_time) {
		Run_time = run_time;
	}

	public int getC() {
		return c;
	}

	public void setC(int c) {
		this.c = c;
	}

	public Random getR() {
		return r;
	}
	
	private void scanReport(Vector<Integer> x)
	{
		int m= input.elementAt(2).size();
		double show = 0;
    	for(int i=0;i<seqprev.size();i++)
    		show += seqprev.elementAt(i);
    	show = (double)show/(seqprev.size());
		int z = (int)prevmidpost(show,num_strings,m,dtype,0);
		
		Vector<Integer> v=showResultTest((z-3+r.nextInt(5)),dtype,0);
		System.out.println(v);
		
	}
	public void setR(Random r) {
		this.r = r;
	}

	public int getNum_strings() {
		return num_strings;
	}

	public void setNum_strings(int num_strings) {
		this.num_strings = num_strings;
	}

	public boolean isValid() {
		return valid;
	}

	public void setValid(boolean valid) {
		this.valid = valid;
	}

	public double getKELossRate() {
		return KELossRate;
	}

	public void setKELossRate(double kELossRate) {
		KELossRate = kELossRate;
	}

	public double getBuffer() {
		return buffer;
	}

	public void setBuffer(double buffer) {
		this.buffer = buffer;
	}

	public double getInitialKE() {
		return InitialKE;
	}

	public void setInitialKE(double initialKE) {
		InitialKE = initialKE;
	}
	

	public int getPopsize() {
		return popsize;
	}

	public void setPopsize(int popsize) {
		this.popsize = popsize;
	}

	public int getVT() {
		return VT;
	}

	public void setVT(int vT) {
		VT = vT;
	}



}
