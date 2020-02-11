package comp559.a2ccd;

import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

/**
 * Particle class for 599 assignment 1
 * @author kry
 */
public class Particle {
    
    /** true means that the particle can not move */
	public boolean pinned = false;
        
    /** The mass of the particle */
	public double mass = 1;
    
    /** current position of the particle */
    public Point2d p = new Point2d();
    
    /** current velocity of the particle */
    public Vector2d v = new Vector2d();
    
    /** initial position of the particle */
    public Point2d p0 = new Point2d();
    
    /** initial velocity of the particle */
    public Vector2d v0 = new Vector2d();
    
    /** force acting on this particle */
    public Vector2d f = new Vector2d();
    
    /**
     * Creates a particle with the given position and velocity
     * @param x
     * @param y
     * @param vx
     * @param vy
     */
    public Particle( double x, double y, double vx, double vy ) {
        p0.set(x,y);
        v0.set(vx,vy);
        reset();
    }
    
    /**
     * Resets the position of this particle
     */
    public void reset() {
        p.set(p0);
        v.set(v0);
        f.set(0,0);
    }
    
    /**
     * Adds the given force to this particle.
     * Note that you probably want to set the force to zero 
     * before accumulating forces. 
     * @param force
     */
    public void addForce( Vector2d force ) {
        f.add(force);
    }
    
    /**
     * Computes the distance of a point to this particle
     * @param x
     * @param y
     * @return the distance
     */
    public double distance( double x, double y ) {
        Point2d tmp = new Point2d( x, y );
        return tmp.distance(p);
    }

    // if the particle si pinned the mass is considered close to infinite, if not it just returns its mass
    public Double getMass(){
        if(pinned==true){
            return Double.POSITIVE_INFINITY;
        }
        return mass;
    }
   
}
