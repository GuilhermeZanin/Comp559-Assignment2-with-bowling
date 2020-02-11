package comp559.a2ccd;

import java.awt.Font;
import java.text.DecimalFormat;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Vector2d;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;

/**
 * Implementation of a simple particle system.
 * 
 * Unlike the first assignment, there is no Function interface for integration.
 * Instead, a symplectic Euler integerator is used.
 * 
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode {
    
    /** the particle list */
    public List<Particle> particles = new LinkedList<Particle>();
    
    /** the spring list (treat this as the edge list for geometry */
    public List<Spring> springs = new LinkedList<Spring>();
    
    /** leaf springs connect 3 particles with a hinge like spring */
    public List<BendingSpring> bendingSprings = new LinkedList<BendingSpring>();
        
    /** The name of the particle system. */
    public String name = "";
    
    boolean grabbed = false;    
    Particle mouseParticle;
    Spring mouseSpring;
    
    public double time = 0;
    
    public RobustCCD robustCCD = new RobustCCD();

    /** temporary variables for stepping */
    private Vector2d tmp = new Vector2d();
    
    /**
     * Creates an empty particle system, and preps the mouse spring
     */
    public ParticleSystem() {
    	mouseParticle = new Particle(0, 0, 0, 0);
    	mouseSpring  = new Spring(mouseParticle, null);
    }
    
    /**
     * Resets the positions of all particles to their initial states
     */
    public void resetParticles() {
        for ( Particle p : particles ) {
            p.reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles, and as such removes all springs too.
     */
    public void clear() {        
        particles.clear();
        springs.clear();
        bendingSprings.clear();
        name = "";
    }    

    /**
     * Updates the velocities with the current forces, given the time step
     * @param h time step
     */
    private void stepVelocities( double h ) {
        for ( Particle p : particles ) {
            if ( p.pinned ) {            
                p.f.set(0,0); // just to make sure!
                p.v.set(0,0);
            } else {
                tmp.scale( h / p.mass, p.f );
                p.v.add( tmp );            
            }
        }
    }

    /**
     * Update the positions with current velocities given the time step 
     * @param h time step
     */
    private void stepPositions( double h ) {
        for ( Particle p : particles ) {
            if ( p.pinned ) continue;
            tmp.scale( h, p.v );
            p.p.add( tmp );
            p.f.set(0,0);
        }
    }
    
    /**
     * Advances time and updates the position of all particles
     * @param h 
     * @return true if update was successful
     */
    public boolean updateParticles( double h ) {

    	// TODO: Observe the main loop here for stepping the system in this assignment!
    	
        computeForces();
        stepVelocities( h );       
        robustCCD.applyRepulsion( h, this );            
        boolean resultOK = robustCCD.check( h, this );
        stepPositions( h );
                        
        time += h;
        return resultOK;
    }
    
    /**
     * Computes forces on all particles
     */
    private void computeForces() {
    	// set the global spring properties
        Spring.ksMod = stiffnessModulation.getValue();
        Spring.ksbMod = springDampingModulation.getValue();
        BendingSpring.kbMod = bendingStiffnessModulation.getValue();
        boolean useCurvature = useCurvatureBendingForce.getValue();
        boolean useStrain = useStrainSpringForce.getValue();
        
        double damping  = viscousDrag.getValue();
        for ( Particle p : particles ) {
            p.f.set( 0, p.mass* (useg.getValue() ? g.getValue() : 0) );
            tmp.scale( -damping*p.mass, p.v );
            p.f.add( tmp );                        
        }
        for ( Spring s : springs ) {
            s.apply( useStrain );
        }
        for ( BendingSpring ls : bendingSprings ) {
            ls.apply( useCurvature );
        }
        if ( grabbed ) {
        	// don't use a strain based compute for mouse spring as
        	// it has zero rest length
        	mouseSpring.ks = mouseSpringStiffness.getValue();
        	mouseSpring.ksb = mouseSpringDamping.getValue();
        	mouseSpring.apply(false);
        }
    }

    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy );
        particles.add( p );
        return p;
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createSpring( Particle p1, Particle p2 ) {
        Spring s = new Spring( p1, p2 ); 
        springs.add( s );         
        return s;
    }
    
    public double[] pack() {
    	double[] state = new double[6*particles.size()];
    	int i = 0;
    	for(Particle p : particles) {
    		state[i++] = p.p.x;
    		state[i++] = p.p.y;
    		state[i++] = p.v.x;
    		state[i++] = p.v.y;
    		state[i++] = p.f.x;
    		state[i++] = p.f.y;
    	}
    	return state;
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    /** Height of the canvas, used by create box in test systems */
    public int height;
    
    /** Width of the canvas, used by create box in test systems */
    public int width;

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        // We'll keep track of the width and the height of the drawable as this 
        // is used to let us easily set up the walls (create box) in the test systems.
        height = drawable.getSurfaceHeight();
        width = drawable.getSurfaceWidth();

        //gl.glColor4d( 0, 0.5, 0.5, 0.75 );
        gl.glColor4d( 0,0,0, 1 );
        gl.glLineWidth( 2 );
        gl.glBegin( GL.GL_LINES );
        for (Spring s : springs) {
            gl.glVertex2d( s.A.p.x, s.A.p.y );
            gl.glVertex2d( s.B.p.x, s.B.p.y );
        }
        gl.glEnd();

        if ( drawParticles.getValue() ) {
            gl.glPointSize( pointSize.getFloatValue() );
            gl.glBegin( GL.GL_POINTS );
            for ( Particle p : particles ) {
                // transparency is used to get smooth edges on the particles
                double alpha = 1;//  0.75;
                if ( p.pinned ) {
                    gl.glColor4d( 1, 0, 0, alpha );
                } else {
                	gl.glColor4d( 0,0,0, alpha );//gl.glColor4d( 0, 0.95,0, alpha );
                }
                gl.glVertex2d( p.p.x, p.p.y );
            }
            gl.glEnd();
        }
    }    
    
    BooleanParameter drawParticles = new BooleanParameter( "draw Particles", true ) ;
    DoubleParameter pointSize = new DoubleParameter("point size", 5, 1, 25);
    
    BooleanParameter useg = new BooleanParameter( "use gravity", true );
    DoubleParameter g = new DoubleParameter( "gravity", 10, 0.01, 100 );
    
    DoubleParameter viscousDrag = new DoubleParameter( "viscous damping", .01, 0, 10 );
    
    BooleanParameter useCurvatureBendingForce = new BooleanParameter( "use curvature bending force", false );
    DoubleParameter bendingStiffnessModulation = new DoubleParameter( "bending stiffness modulation", 1e3, 100, 1e5 );
    
    BooleanParameter useStrainSpringForce = new BooleanParameter( "use strain spring force", false );
    DoubleParameter stiffnessModulation = new DoubleParameter( "spring stiffness modulation", 100, 0.01, 100000 );
    DoubleParameter springDampingModulation = new DoubleParameter( "spring damping modulation", 1, 0, 10 );
    
    DoubleParameter mouseSpringStiffness = new DoubleParameter( "mouse spring stiffness", 0.3, 0.1, 10 );
    DoubleParameter mouseSpringDamping = new DoubleParameter( "mouse spring stiffness", 10, 0.1, 100 );
    
                
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
                
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder( new TitledBorder("Simulation Parameters"));
        ((TitledBorder) vfp1.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp1.add( useg.getControls() );
        vfp1.add( g.getSliderControls(true) );
        vfp1.add( viscousDrag.getSliderControls(false) );
        vfp1.add( useCurvatureBendingForce.getControls() );
        vfp1.add( bendingStiffnessModulation.getSliderControls(true) );
        vfp1.add( useStrainSpringForce.getControls() );
        vfp1.add( stiffnessModulation.getSliderControls(true) );
        vfp1.add( springDampingModulation.getSliderControls(false) );
        vfp1.add( mouseSpringStiffness.getSliderControls(true));
        vfp.add( vfp1.getPanel() );
        
        vfp.add( robustCCD.getControls());

        return vfp.getPanel();        
    }

    private DecimalFormat df = new DecimalFormat("0.000");

    @Override
    public String toString() {
        String s = name + "\n" +
 	           "CCD iters = " + robustCCD.iters + "\n" +
        	   "particles = " + particles.size() + "\n" + 
 	           "time = " + df.format(time) + "\n" +
//	           "stiffness = " + stiffnessModulation.getValue() + "\n" +               
//	           "spring damping = " + springDampingModulation.getValue() + "\n" +
//	           "viscous damping = " + viscousDrag.getValue() + "\n" +
//	           "bending stiffness = " + bendingStiffnessModulation.getValue() + "\n" +
//	           "restitution = " + robustCCD.restitutionValue.getValue() + "\n"+ 
//	           "H = " + robustCCD.minDist.getValue() + "\n" +
	           (robustCCD.useJacobi.getValue() ? "using Jacobi\n" : "");
        return s;
    }
    
}
