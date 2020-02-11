package comp559.a2ccd;

import javax.vecmath.Vector2d;

/**
 * Simple spring class where spring forces are based on strain
 * @author kry
 */
public class Spring {

    Particle A = null;
    
    Particle B = null;
    
    /** spring stiffness modulation */
    static double ksMod = 1;

    /** spring stiffness modulation */
    static double ksbMod = 1;

    /** Spring stiffness */
    double ks = 1;
    
    /** spring damping */
    double ksb = 1;
    
    /** Rest length */
    double l0 = 0;
    
    /**
     * Creates a spring connecting two particles.
     * The rest length should be set
     * @param p1
     * @param p2
     */
    public Spring( Particle p0, Particle p1 ) {
        this.A = p0;
        this.B = p1;
        if ( p0 == null || p1 == null ) return;
        computeRestLength();
    }
    
    /**
     * Computes the rest length of the connected particles
     */
    public void computeRestLength() {
        l0 = A.p0.distance( B.p0 );
    }
    
    /**
     * Sets the rest length of the connected particles with their current positions
     */
    public void setRestLength() {
        l0 = A.p.distance( B.p );
    }
    
    public void apply( boolean strainBased ) {
    	if ( strainBased ) {
    		applyStrainBasedSpring();
    	} else {
    		applyLengthBasedSpring();
    	}
    }
    
    /**
     * Applies spring length based force and damping forces to the two particles
     */
    public void applyLengthBasedSpring() {
        
        force.sub( B.p, A.p );
        double length  = force.length();
        force.normalize();
        force.scale( (length-l0) * ks * ksMod );
        A.addForce(force);
        force.scale(-1);
        B.addForce(force);
        
        force.sub( B.p, A.p );
        force.normalize();
        v.sub(B.v, A.v);
        double rv = force.dot(v);
        force.scale( rv * ksb * ksbMod );
        A.addForce(force);
        force.scale(-1);
        B.addForce(force);            
    }

    
    private Vector2d force = new Vector2d();
    private Vector2d v = new Vector2d();

    /**
     * Applies spring strain based force and damping forces to the two particles
     */
    public void applyStrainBasedSpring() {
        
        force.sub( B.p, A.p );
        double strain  = force.length() / l0 - 1;
        force.normalize();
        force.scale( strain * ks * ksMod );
        A.addForce(force);
        force.scale(-1);
        B.addForce(force);
        
        force.sub( B.p, A.p );
        force.normalize();
        v.sub(B.v, A.v);
        double rv = force.dot(v);
        force.scale( rv * ksb * ksbMod );
        A.addForce(force);
        force.scale(-1);
        B.addForce(force);            
    }
    
}
