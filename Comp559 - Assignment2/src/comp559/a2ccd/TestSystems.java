package comp559.a2ccd;

import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;

/**
 * Somewhat ugly class to procedurally help in creating a variety of interesting test scenes
 * 
 * @author kry
 */
public class TestSystems {

    public BooleanParameter runAlphabetFactory = new BooleanParameter( "run alphabet soup factory!", false );
    
    public BooleanParameter runPastaFactory = new BooleanParameter( "run pasta factory!", false );
    
    /**
     * Anything less than 100 is likely to end up with letters being generated over top one another and killing 
     * the simulation
     */
    private IntParameter interval = new IntParameter( "new letter step interval", 130, 100, 200 );

    private IntParameter stringParticleInterval = new IntParameter( "new pasta particle step interval", 3, 1, 40 );
        
    private DoubleParameter initialVelocity = new DoubleParameter( "initial velocity", 50, -75, 100 );
    
    private BooleanParameter clearFirst = new BooleanParameter( "clear current system before creating new systems", true );
    
    private AlphabetSoupFactory alphabetSoupFactory;

    private ParticleSystem system;
    
    private DoubleParameter mass = new DoubleParameter("p4p5 mass for sandwitch", 1, 1, 1000);
    
    // for convenience we'll keep a copy of the particles, springs, and leaf springs inside a system,
    // though this is a bit gross    
    
    private List<Particle> particles;
    private List<Spring> springs;
    private List<BendingSpring> bendingSprings;
    
    /**
     * Creates a new test system 
     * @param system
     */
    public TestSystems( ParticleSystem system ) {
        this.system = system;
        particles = system.particles;
        springs = system.springs;
        bendingSprings = system.bendingSprings;
        alphabetSoupFactory = new AlphabetSoupFactory( system );
    }
    
    /**
     * Quick and dirty generic test generation button
     * @author kry
     */
    private class TestButton extends JButton implements ActionListener {
        private static final long serialVersionUID = 1L;
        private int testNumber;
        public TestButton( String name, int testNumber ) {
            super( name );
            this.testNumber = testNumber;
            addActionListener( this );
        }
        @Override
        public void actionPerformed(ActionEvent e) {
            createSystem(this.testNumber);
        }
    }
    
    /**
     * Gets the control panel for setting different systems.
     * @return the control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();

        JButton addBox = new JButton("Add a open box" );
        addBox.addActionListener( new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
	            createBox();
			}
		});
        vfp.add( addBox );
        vfp.add( clearFirst.getControls() );
        
        
        JPanel buttons = new JPanel();
        buttons.setLayout(new GridBagLayout());
        for ( int i = 0; i < tests.length; i++ ) {   
        	JButton button = new TestButton(tests[i], i);
        	GridBagConstraints c = new GridBagConstraints();
        	c.fill = GridBagConstraints.HORIZONTAL;
        	c.weightx = 0.33;
        	c.gridx = i%3;
        	c.gridy = i/3;
        	buttons.add(button, c);
        }

        vfp.add( buttons );

        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder( new TitledBorder("Endless Pasta Factory controls"));
        ((TitledBorder) vfp1.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp1.add( runPastaFactory.getControls() );                
        vfp1.add( stringParticleInterval.getControls() );
        vfp1.add( initialVelocity.getSliderControls(false) );
        vfp.add( vfp1.getPanel() );

        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Alphabet Soup Factory Controls"));
        ((TitledBorder) vfp0.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp0.add( runAlphabetFactory.getControls() );
        vfp0.add( interval.getControls() );
        vfp.add( vfp0.getPanel() );

        vfp.add( alphabetSoupFactory.getControls() );

        VerticalFlowPanel vfp2 = new VerticalFlowPanel();
        vfp2.setBorder( new TitledBorder("Sandwich Momentum Transfer Test Controls"));
        ((TitledBorder) vfp0.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp2.add( mass.getSliderControls(false) );
        vfp2.add( new TestButton( "Sandwich", 100 ) ); // special ID, not included in the named list!
        vfp.add( vfp2.getPanel() );
        
        return vfp.getPanel();   
    }
    
    /** 
     * Creates a pinned box with an open top that matches the current display window
     */
    private void createBox() {
        double h = system.height;
        double w = system.width;
        Particle p1 = system.createParticle(     5,     5, 0, 0 ); p1.pinned = true;
        Particle p2 = system.createParticle(     5, h - 5, 0, 0 ); p2.pinned = true;
        Particle p3 = system.createParticle( w - 5, h - 5, 0, 0 ); p3.pinned = true;
        Particle p4 = system.createParticle( w - 5,     5, 0, 0 ); p4.pinned = true;        
        system.createSpring(p1, p2);
        system.createSpring(p2, p3);
        system.createSpring(p3, p4);
    }

    /**
     * Creates a free particle above a pinned line segment
     */
    private void createSimple() {
	    Particle p1 = new Particle( 100, 400, 0, 0 );
		Particle p2 = new Particle( 500, 400, 0, 0 );
	    particles.add( p1 );
	    particles.add( p2 );
	    p1.pinned = true;
	    p2.pinned = true;
	    springs.add( new Spring( p1, p2 ) );
		Particle p3 = new Particle( 300, 50, 0, 0 );
		particles.add(p3);
    }
    
    private void createPendulumPair() {
    	Particle p1 = new Particle( 320, 100, 0, 0 );
        Particle p2 = new Particle( 520, 100, 0, 0 );
        particles.add( p1 );
        particles.add( p2 );
        p1.pinned = true;
        springs.add( new Spring( p1, p2 ) );
        Particle p3 = new Particle( 300, 150, 0, 0 );
        Particle p4 = new Particle( 260, 150, 0, 0 );
        particles.add( p3 );
        particles.add( p4 );
        p3.pinned = true;
        springs.add( new Spring( p3, p4 ) );        
    }

    /**
     * Creates a vertical chain, haning from the specified location, with given nodes and length.
     * Masses and set according to spring lengths.
     * Does not include bend springs, but perhaps it should!
     * @param xpos
     * @param ypos
     * @param N
     * @param L
     */
    private void createChain( double xpos, double ypos, int N, double L ) {
    	double segLength = L/N;
    	Particle p1, p2;
        p1 = new Particle( xpos, ypos, 0, 0 );
        p1.mass = 0;
        p1.pinned = true;
        particles.add( p1 );
        for ( int i = 0; i <= N; i++ ) {
            ypos += i==0?20:segLength;
            p2 = new Particle( xpos, ypos, 0, 0 );
            p2.mass = 0;
            particles.add( p2 );
            Spring s = new Spring( p1, p2 );
            springs.add( s );
            p1.mass += s.l0/20/2;
            p2.mass += s.l0/20/2;
            p1 = p2;
        }
    }
    
    /**
     * Creates a row of chains with different numbers of elements.
     * Would desire the same behaviour!
     */
    private void createChains() {
    	int n = 5;
    	for ( int i = 0; i < 10; i ++ ) {
            createChain( 320 + i*10, 100, n++, 200 );
    	}
    }

    private void createTriangleTruss() {
	    Point2d p = new Point2d(100, 100);
	    Vector2d d = new Vector2d(20, 0);
	    Particle p1, p2, p3, p4;
	    p1 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
	    particles.add(p1);
	    p2 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
	    particles.add(p2);
	    springs.add(new Spring(p1, p2));
	    p1.pinned = true;
	    p2.pinned = true;
	    p.add(d);
	    p.add(d);
	    int N = 10;
	    for (int i = 1; i < N; i++) {
	        // d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );
	        p3 = new Particle(p.x - d.y, p.y + d.x, 0, 0);
	        p4 = new Particle(p.x + d.y, p.y - d.x, 0, 0);
	        particles.add(p3);
	        particles.add(p4);
	        springs.add(new Spring(p3, p1));
	        //springs.add(new Spring(p3, p2));
	        springs.add(new Spring(p4, p1));
	        springs.add(new Spring(p4, p2));
	        springs.add(new Spring(p4, p3));
	        p1 = p3;
	        p2 = p4;
	
	        p.add(d);
	        p.add(d);
	    }
	}
    
    /**
     * Creates a zig zag of particles with bending springs in between,
     * a good test of non zero rest angle or curvature.
     */
    private void createZigZag() {
    	 int N = 20;
         int xpos = 100;
         Particle p0, p1, p2;
         
         p0 = null;
         p1 = null;
         p2 = null;            
         for ( int i = 0; i < N; i++ ) {               
             p2 = new Particle( xpos, 100 + 20*(i%2), 0, 0 );                
             particles.add( p2 );
             if ( i < 2 ) p2.pinned = true;
             if ( p1 != null ) springs.add( new Spring( p1, p2 ) );
             if ( p0 != null ) bendingSprings.add( new BendingSpring( p0, p1, p2 ) );
             p0 = p1;
             p1 = p2;                
             xpos += 20;
         }
    }
    
    /**
     * Like a chain, but creates with bending springs too!
     * @param xpos
     * @param ypos
     * @param dx
     * @param dy
     */
    private void createHair( double xpos, double ypos, double dx, double dy, int N ) {
    	Particle p0 = null;
    	Particle p1 = null;
        Particle p2 = null;            
        for ( int i = 0; i < N; i++ ) {               
            p2 = new Particle( xpos + dx*i, ypos+ dy*i, 0, 0 );                
            particles.add( p2 );
            if ( i < 2 ) p2.pinned = true;
            if ( p1 != null ) springs.add( new Spring( p1, p2 ) );
            if ( p0 != null ) bendingSprings.add( new BendingSpring( p0, p1, p2 ) );
            p0 = p1;
            p1 = p2;                
        }
    }
    
    private void createHairsHard() {
    	int N = 20;
        int M = 5;
        int offset = 10;
        for ( int j = 0; j < M; j++ ) {
        	createHair( 100              , 100 + offset*2*j,           20, 0, 20);
        	createHair( 100 + (N-2) * 20 , 100 + offset*2*j + offset, -20, 0, 20);
        }
    }
    
    private void createHairsEasy() {
    	int N = 20;
        int M = 5;
        int offset = 15;
        for ( int j = 0; j < M; j++ ) {
        	createHair( 100              , 100 + offset*2*j,           20, 0, 20);
        	createHair( 100 + (N-2) * 20 , 100 + offset*2*j + offset, -20, 0, 20);
        }
    }
    
    private void createBulletAndLayers() {
    	int numHairs = 10;
    	double spacing = 8;
    	int N = 15;
    	for ( int i = 0; i < numHairs; i++ ) {
    		createHair( 300+i*spacing, 100, 0, 10, N);
    	}
    	Particle p = new Particle( 200, 150, 20, 0 );
    	p.mass = 300;
    	particles.add( p );
    	system.useg.setValue(false);
    	system.bendingStiffnessModulation.setValue( 1e3 );
    	// should change the base stiffness of the hairs instead!
    }
    
    /** 
     * Creates a test in which only symmetry preserved by Jacobi will lead to a 
     * balanced result.
     */
    private void createJacobiTest() {
    	Particle p1, p2;
        // left end
        p1 = new Particle(250, 400, 0, 0 );
        p1.pinned = true;
        particles.add( p1 );
        p2 = new Particle(250, 250, 0, 0 );
        particles.add( p2 );
        springs.add( new Spring(p1, p2) );
        
        // right end
        p1 = new Particle(500, 400, 0, 0 );
        p1.pinned = true;
        particles.add( p1 );
        p2 = new Particle(500, 250, 0, 0 );
        particles.add( p2 );
        springs.add( new Spring(p1, p2) );
        
        // horizontal
        p1 = new Particle(250, 100, 0, 0 );
        particles.add( p1 );
        p2 = new Particle(500, 100, 0, 0 );
        particles.add( p2 );
        springs.add( new Spring(p1, p2) );  
    }
    
    public void createElasticRodCollision() {
    	// create a rod and set its velocity
    	int N = 20;
    	int xpos = 30;
    	int ypos = 300;
    	float dx = 100f/N;
        Particle p1 = new Particle(xpos, ypos, 10, 0 );
    	particles.add( p1 );
        xpos += dx;
    	for ( int i = 0; i < N; i++ ) {
            Particle p2 = new Particle(xpos, ypos, 10, 0 );
            particles.add( p2 );
            springs.add( new Spring(p1, p2) );
            p1 = p2;
            xpos += dx;
        }
    	
    	xpos = 300;
    	p1 = new Particle(xpos, ypos, 0, 0 );
    	particles.add( p1 );
        xpos += dx;
    	for ( int i = 0; i < N; i++ ) {
            Particle p2 = new Particle(xpos, ypos, 0, 0 );
            particles.add( p2 );
            springs.add( new Spring(p1, p2) );
            p1 = p2;
            xpos += dx;
        }
    	system.useg.setValue(false);
    }
    
    private void createSimpleBendTest() {
	 	int xpos = 300;
		int ypos = 300;
		float dx = 100;
	    Particle p0 = new Particle(xpos, ypos, 0, 0 );
	    p0.pinned = true;
	    Particle p1 = new Particle(xpos+dx, ypos, 0, 0 );
	    p1.pinned = true;
	    Particle p2 = new Particle(xpos, ypos+dx, 0, 0 );
	    particles.add( p0 );
	    particles.add( p1 );
	    particles.add( p2 );
	    Spring s1 = new Spring( p0, p1 );
	    springs.add(s1);
	    Spring s2 = new Spring( p1, p2 );
	    springs.add(s2);
	    BendingSpring bs = new BendingSpring(p0, p1, p2);
	    system.bendingSprings.add(bs);
	}

    /**
     * Creates a particle with velocity between two 
     * springs, with one pinned, and the other pair being potentially
     * quite massive.
     * How many iterations of resolution will be necessary?
     * Note the coefficient of restitution is set to 1 here!
     */
    private void createSandwich() {
    	Particle p1 = new Particle( 100, 200, 0, 0 );
    	Particle p2 = new Particle( 100, 400, 0, 0 );
    	Particle p3 = new Particle( 100.1, 300, 200, 0 );
    	Particle p4 = new Particle( 100.2, 200, 0, 0 );
    	Particle p5 = new Particle( 100.2, 400, 0, 0 );
        particles.add( p1 );
        particles.add( p2 );
        particles.add( p3 );
        particles.add( p4 );
        particles.add( p5 );            
        p1.pinned = true;
        p2.pinned = true;
        p3.mass = 2;
        p4.mass = mass.getValue();
        p5.mass = mass.getValue();
        springs.add( new Spring( p1, p2 ) );
        springs.add( new Spring( p4, p5 ) );
        system.robustCCD.restitutionValue.setValue(1);
        system.robustCCD.repulsion.setValue(false);
        system.useg.setValue(false);
    }
    /**
     * 
     * @param cx
     * @param cy
     * @param r
     * @param N
     * @param k  create k of N points... close if k = N
     * @param thetaOffset if zero starts at top
     * @param pinned if true will do the first two
     */
    private void createCircle( double cx, double cy, double r, int N, int k, double thetaOffset, boolean pinned ) {
    	Particle p0=null, p1=null, p2=null, first=null, second=null;
    	for ( int i = 0; i < k; i++ ) {
    		double theta = i*2*Math.PI/N;
    		p2 = new Particle( cx - Math.sin(theta + thetaOffset )*r, cy - Math.cos(theta + thetaOffset )*r, 0, 0 );  
    		p2.mass = 0; 
     		if ( i == 0 ) first = p2;
     		if ( i == 1 ) second = p2;
    		particles.add( p2 );
    		if ( i<2 && pinned ) p2.pinned = true;
    		if ( p1 != null ) {
    			Spring s = new Spring( p1, p2 );
    			springs.add( s );
    	    	p1.mass += s.l0/20/2;
    	        p2.mass += s.l0/20/2;        

    		}
            if ( p0 != null ) bendingSprings.add( new BendingSpring( p0, p1, p2 ) );
            p0 = p1;
            p1 = p2;        
    	}
    	if ( k == N ) {
    	Spring s = new Spring (p1, first); 
    	springs.add( s );
    	p1.mass += s.l0/20/2;
        first.mass += s.l0/20/2;        
    	bendingSprings.add( new BendingSpring(p0, p1, first ) );
    	bendingSprings.add( new BendingSpring(p1, first, second) );
    	}
    }

    private void createCircles() {
    	createCircle(150, 150, 50, 20, 20, 0, false);
    	createCircle(300, 150, 50, 30, 30, 0, false);
    	createCircle(450, 150, 50, 40, 40, 0, false);
    	createBox();
    }
    
    /** How about a hairy ball ? */
    private void createBentHairs() {
    	int N = 20;
    	double cx = 300;
    	double cy = 200;
    	double r = 50;
    	Particle[] circlePoints = new Particle[N];
    	for ( int i = 0; i < N; i++ ) {
     		double theta = i*2*Math.PI/N;
     		circlePoints[i] = new Particle( cx - Math.sin(theta)*r, cy - Math.cos(theta)*r, 0, 0 );                
     		particles.add( circlePoints[i] );
    	}
    	
    	int k = 10;
    	for ( int i = 0; i < N; i++ ) {
     		springs.add( new Spring( circlePoints[i], circlePoints[(i+1)%N] ) );
     		bendingSprings.add( new BendingSpring( circlePoints[i], circlePoints[(i+1)%N], circlePoints[(i+2)%N] ) );
     		double theta = i*2*Math.PI/N;     			
     		Particle p0 = circlePoints[(i+N-1)%N];
     		Particle p1 = circlePoints[i];
     		for ( int j = 1; j < k; j++ ) {
         		double px = cx - Math.sin(theta)*(r+j*12);
 				double py = cy - Math.cos(theta)*(r+j*12);
 				Particle p2 = new Particle( px, py, 0, 0 );
 				particles.add( p2 );
 				springs.add( new Spring( p1, p2 ) );
 				bendingSprings.add( new BendingSpring( p0, p1, p2 ) );
 	     		p0 = p1;
 	     		p1 = p2;
     		}            
     	}
     	createBox();
     }
    
    /** Bowling! */
    private void createBowling() {
    	// ground
    	Particle p1 = new Particle( 0, 450, 0, 0 );
		Particle p2 = new Particle( 640, 450, 0, 0 );
	    particles.add( p1 );
	    particles.add( p2 );
	    p1.pinned = true;
	    p2.pinned = true;
	    springs.add( new Spring( p1, p2 ) );
	    
	    //pins
	    //Adds 3 pins in a do while loop
	    int i=0;
	    int moveX=1;
	    do{
	    	i++;
		//  which particle to which position
		//		    11 12
		//		   9   10
		//		  7     8
		//		  5     6
		//		    3  4
	    	Particle p3 = new Particle( 20+moveX, 440, 0, 0 );
			Particle p4 = new Particle( 30+moveX, 440, 0, 0 );
			Particle p5 = new Particle( 10+moveX, 430, 0, 0 );
			Particle p6 = new Particle( 40+moveX, 430, 0, 0 );
			Particle p7 = new Particle( 5+moveX, 410, 0, 0 );
			Particle p8 = new Particle( 45+moveX, 410, 0, 0 );
			Particle p9 = new Particle( 5+moveX, 395, 0, 0 );
			Particle p10 = new Particle( 45+moveX, 395, 0, 0 );
			Particle p11 = new Particle( 20+moveX, 365, 0, 0 );
			Particle p12 = new Particle( 30+moveX, 365, 0, 0 );
			moveX+=60;
			
			
			//springs to keep the pin together
			//bottom, left side and top
		    springs.add( new Spring( p3, p4 ) );
		    springs.add( new Spring( p3, p5 ) );
		    springs.add( new Spring( p5, p7 ) );
		    springs.add( new Spring( p7, p9 ) );
		    springs.add( new Spring( p9, p11 ) );
		    springs.add( new Spring( p11, p12 ) );
		    //right side
		    springs.add( new Spring( p4, p6 ) );
		    springs.add( new Spring( p6, p8 ) );
		    springs.add( new Spring( p8, p10 ) );
		    springs.add( new Spring( p10, p12 ) );
			
			particles.add( p3 );
		    particles.add( p4 );
		    particles.add( p5 );
		    particles.add( p6 );
		    particles.add( p7 );
		    particles.add( p8 );
		    particles.add( p9 );
		    particles.add( p10 );
		    particles.add( p11 );
		    particles.add( p12 );
	    }while(i<=2);
	    
	    //Ball - being a poligon with 10 sides with a center point
	    Particle pTL = new Particle( 440, 380, -75, 0 ); // top left
	    Particle pTR = new Particle( 480, 380, -75, 0 ); // top right
	    Particle pML = new Particle( 430, 400, -75, 0 ); // middle left
	    Particle pMR = new Particle( 490, 400, -75, 0 ); // middle right
	    Particle pCT = new Particle( 460, 370, -75, 0 ); // center Top
	    Particle pC = new Particle( 460, 400, -75, 0 ); // center
	    Particle pCB = new Particle( 460, 430, -75, 0 ); // center bottom
	    Particle pBL = new Particle( 440, 420, -75, 0 ); // bottom left
	    Particle pBR = new Particle( 480, 420, -75, 0 ); // bottom right
	    
	    //linking them clockwise the borders
	    springs.add( new Spring( pTL, pCT ) );
	    springs.add( new Spring( pCT, pTR ) );
	    springs.add( new Spring( pTR, pMR ) );
	    springs.add( new Spring( pMR, pBR ) );
	    springs.add( new Spring( pBR, pCB ) );
	    springs.add( new Spring( pCB, pBL ) );
	    springs.add( new Spring( pBL, pML ) );
	    springs.add( new Spring( pML, pTL ) );
	    
	    
	    particles.add( pTL );
	    particles.add( pTR );
	    particles.add( pML );
	    particles.add( pMR );
	    particles.add( pCT );
	    particles.add( pC );
	    particles.add( pCB );
	    particles.add( pBL );
	    particles.add( pBR );
	    
	    system.useg.setValue(false);
    }
    
    public String[] tests = {            
            "bullet+layers",
            "simple",
            "chains",
            "pendulum pair",
            "zig zag",
            "interwoven hairs hard",
            "interwoven hairs easy",
            "triangle truss",
            "vertical chain", 
            "Jacobi test",
            "elastic rods",
            "bend test",
            "Circles",
            "bent hairs",
            "Bowling ball!",
            "<empty>",
            "<empty>",
            "<empty>",
    };
    
    /**
     * Creates one of a number of simple test systems.
     *
     * Small systems are more useful for debugging!
     * 
     * @param which
     */
    private void createSystem( int which ) {
        if ( clearFirst.getValue() ) {
        	system.clear();
        }       
        system.name = tests[which];
        
        system.useg.setValue(true); // wanted for most of these examples
        system.robustCCD.restitutionValue.setValue(0);
        system.robustCCD.repulsion.setValue(true);

        if (which == 0 ) {
        	createBulletAndLayers();
        } else if ( which == 1 ) {
        	createSimple();
    	} else if ( which == 2) { 
            createChains();
        } else if ( which == 3) {
            createPendulumPair();
        } else if ( which == 4 ) {            
        	createZigZag(); 
        } else if ( which == 5 ) {
        	createHairsHard();           
        } else if ( which == 6 ) {
        	createHairsEasy();
        } else if ( which == 7 ) {
        	createTriangleTruss();
        } else if ( which == 8 ) {
            createChain( 320, 100, 10, 200 );
        } else if ( which == 9 ) {   
        	createJacobiTest();
        } else if ( which == 10 ){
        	createElasticRodCollision();
        } else if ( which == 11 ) {
        	createSimpleBendTest();        
        } else if ( which == 12 ) { 
        	createCircles();
        } else if ( which == 13 ) {
        	createBentHairs();
        } else if ( which == 14 ) {
        	// TODO: create your own test systems here!
        	// set the names for the interface buttons in the String array "tests" above
        	createBowling();
        } else if ( which == 15 ) {

        } else if ( which == 16 ) {

        } else if ( which == 17 ) {

        } else if ( which == 100 ) { 
        	createSandwich();
        }        
    }
    
    /**
     * reset the factory so that it always starts with the first letter
     */
    public void resetFactory() {
        stepCount = 0;
        letterIndex = 0;
        prevParticle0 = null;
        prevParticle1 = null;
    }
    
    Particle prevParticle0 = null;
    Particle prevParticle1 = null;
    
    private int stepCount = 0;
    
    public final String alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
    public int letterIndex = 0;
    
    /**
     * Generates letters at regular intervals if the factory is turned on.
     */
    public void step() {
        if ( ! runAlphabetFactory.getValue() && ! runPastaFactory.getValue() ) return;
            
        stepCount++;
        
        if ( runAlphabetFactory.getValue() ) {
            if (stepCount % interval.getValue() == 1) {
                alphabetSoupFactory.createLetter(system, "" + alphabet.charAt(letterIndex++%alphabet.length()), 50, -20 );
            }
        }
        
        if ( runPastaFactory.getValue() ) {
            if ( stepCount % stringParticleInterval.getValue() == 1 ) {
                Particle p = new Particle( 100, 30, initialVelocity.getValue(), Math.sin(stepCount/20.0)* initialVelocity.getValue() );
                particles.add( p );
                if ( prevParticle1 != null ) {
                    Spring s = new Spring( prevParticle1, p );
                    s.setRestLength(); // set rest length with current positions.
                    springs.add( s );
                }
                if ( prevParticle0 != null && prevParticle1 != null ) {
                    BendingSpring bs = new BendingSpring(prevParticle0,prevParticle1,p);
                    bs.theta0 = 0;
                    bs.kappa0 = 0;
                    bendingSprings.add( bs );
                }
                prevParticle0 = prevParticle1;
                prevParticle1 = p;
            }
        }
        
    }
    
}
