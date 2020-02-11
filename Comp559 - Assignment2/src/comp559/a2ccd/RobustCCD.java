package comp559.a2ccd;

import java.awt.Font;
import java.util.Vector;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.Vector2d;
import javax.xml.bind.util.ValidationEventCollector;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import org.netlib.lapack.Dlaed6;

/**
 * Implementation of a robust collision detection
 * @author kry
 */
public class RobustCCD {
	        
    /** number of iterations in the last CCD processing loop, to keep an eye on how tricky the problem is */
    int iters;
    
    /**
     * Creates the new continuous collision detection and response object
     */
    public RobustCCD() {
        // do nothing
    }
    
    /** Might want to turn off collisions when testing? */
    BooleanParameter collision = new BooleanParameter( "apply collision impulses", true );

    /** Ignore this parameter unless you want to explore a Jacobi type resolution of collision */
    BooleanParameter useJacobi = new BooleanParameter( "use Jacobi", false );

    /** Use this as the maximum number of iterations, feel free to modify default, or the maximum! */
    IntParameter maxIterations = new IntParameter("maximum iterations", 60, 30, 300 );

    BooleanParameter repulsion = new BooleanParameter( "apply repulsion impulses", true );

    DoubleParameter restitutionValue = new DoubleParameter( "restitution", .0001, 0, 1 );

    DoubleParameter minDist = new DoubleParameter( "min distance (H)", 2, 0.1, 10 );

    public JPanel getControls() {
    	VerticalFlowPanel vfp = new VerticalFlowPanel();
    	vfp.setBorder( new TitledBorder("Robust CCD Controls"));
        ((TitledBorder) vfp.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
    	vfp.add( collision.getControls() );
    	vfp.add( useJacobi.getControls() );
    	vfp.add(maxIterations.getSliderControls());
    	vfp.add( repulsion.getControls() );
    	vfp.add( restitutionValue.getSliderControls(false) );
    	vfp.add( minDist.getSliderControls(true));    	
    	return vfp.getPanel();
    }
    
    /**
     * Try to deal with contacts before they happen
     * @param h
     * @param system
     */
    public void applyRepulsion( double h, ParticleSystem system ) {
    	if ( ! repulsion.getValue() ) return;
    	// TODO: apply repulsion forces
    	// use minDist.getValue() as the thickness
    	// use your spring stiffness for the repulsion stifness, or create 
    	// new parameters and set their value!
    	
    	//doing the repulsion forces (spring based repulsion forces)


        //Do for every particle in the partyclesystem
        for(Particle tempParticle : system.particles) {

            for (Spring tempSpring : system.springs) {

                if(tempParticle == tempSpring.A || tempParticle == tempSpring.B)
                    continue;

                Vector2d normal = computeNormal(tempSpring.A, tempSpring.B, tempParticle, h);
                Vector2d distanceCB = new Vector2d();

                //getting vector distance
                distanceCB.sub(tempParticle.p, tempSpring.B.p);
                distanceCB.scaleAdd(h, tempParticle.v, distanceCB);
                distanceCB.scaleAdd(-h, tempSpring.B.v, distanceCB);

                //d
                double distance=minDist.getValue()-distanceCB.dot(normal); // gives a value, should be a minus but the disntance dotted with normal is negative, so a plus

                //if distance is less or euqal to 0 its on the space set
                if(distance<=0)
                    continue;

                Double alpha = computeAlpha(tempSpring.A, tempSpring.B, tempParticle, h);

                //computer L
                double l = distanceCB.length();
                //checks distance of the particles around, particle particle repulsion
                if(Double.isFinite(alpha) && alpha*l>=-minDist.getValue() && alpha*l<=l+minDist.getValue()){
                        Vector2d relVelocity = computeRelVelocity(tempSpring.A, tempSpring.B, tempParticle, alpha);

                        // paper Ir equation / temporary equation 0.1d/deltat - vn
                        Double tempVariableInside=tempParticle.getMass()*(((0.1*distance)/h)+relVelocity.dot(normal));

                        if(tempVariableInside<=0)
                            continue;

                        Double repulsionImpulse = Math.min(h*tempSpring.ks*distance, tempVariableInside);

                        distributeImpulse(tempSpring.A, tempSpring.B,tempParticle, repulsionImpulse, normal, alpha);

                }

            }
        }
    	

    }
    
    /**
     * Checks all collisions in interval t to t+h
     * @param h
     * @param system 
     * @return true if all collisions resolved
     */
    public boolean check( double h, ParticleSystem system ) {        
    	if ( ! collision.getValue() ) return true; // pretend everything is OK!
    	
    	// TODO: find collisions and apply impulses to resolve them
    	// use maxIterations.getValue() as max iteraitons before giving up
    	// use restitutionValue.getValue() for computing the impulse

    	//creates a boolean to check if there is collision
        boolean collision = false;

        //iterations went thro the loop
        iters=0;
        // add epsilon to h to make it more robust
        double epsilon = 1/1000;
        
        //enters the loop to check for colisions and adding impulses 
        //and continues while colisions is true or the number of iterations is less than the maximum iterations set
        do{
            //increasses the number of iteratiosn and exits the loop if there is colision or the iterations exceeded the max iteration value
            iters++;

            collision = false;
            //Do for every particle in the partyclesystem
            for(Particle tempParticle : system.particles) {

                for(Spring tempSpring : system.springs) {

                    // if the particle colliding is the same as particle A and B on spring, skips loop and goes to the next step
                    if(tempParticle == tempSpring.A || tempParticle == tempSpring.B)
                        continue;

                    // if Particle C is not particle A and B, calculate the time, tempParticle is particle C
                    Double t  = computeT(tempSpring.A, tempSpring.B,tempParticle);
                    if(Double.isFinite(t) && t>0 && t<=h+epsilon){// chekc if T is not infinite and bigger than 0, no colisions in the past, and has to be smaller than the timestep (h)
                        double alpha = computeAlpha(tempSpring.A, tempSpring.B, tempParticle, t);
                        // check if alpha is finite, bigger than 0, and smaller/qual to 1+epsilon // makes sure it doesnt skip particles (doesnt hit either A or B)
                        if(Double.isFinite(alpha) && alpha>0-epsilon && alpha<=1+epsilon) {
                            // calculate the normal
                            Vector2d normal = computeNormal(tempSpring.A, tempSpring.B, tempParticle, t);

                            // calculates impulse
                            Double impulse = computeImpulse(tempSpring.A, tempSpring.B, tempParticle, alpha, normal);

                            // adds impulse to all the particles
                            distributeImpulse(tempSpring.A, tempSpring.B, tempParticle, impulse, normal, alpha);
                            // adds colision
                            collision=true;
                        }
                    }
                }

            }
        }while(collision == true && iters<maxIterations.getValue());


        
        return !collision;
    }

    // calculated t, when there is colision, when the 3 particles are in a line/colinear
    // determinant of |C-A|
    //                |B-A| = 0
    public double computeT(Particle A, Particle B, Particle C){
        Vector2d alpha = new Vector2d(); //position of colision
        Vector2d beta = new Vector2d();
        Vector2d alphaDot = new Vector2d(); // velocites
        Vector2d betaDot = new Vector2d();

        alpha.sub(B.p, A.p); // sub of positions, sub does not matter order, because are determinants
        beta.sub(B.p, C.p); // sub of positions

        alphaDot.sub(B.v, A.v); // sub of velocities, sub does not matter order, because are determinants
        betaDot.sub(B.v, C.v); // sub of velocities

        double a, b, c; // finding the ts fo the quadratic formula a+b+c=0 / t2+t1+t0 = 0

        a=alphaDot.x*betaDot.y-alphaDot.y*betaDot.x;
        b=alpha.x*betaDot.y+alphaDot.x*beta.y-alpha.y*betaDot.x-alphaDot.y*beta.x;
        c=alpha.x*beta.y-alpha.y*beta.x;

        double t1, t2;
        // if a is zeor then the t is going to be t=-b/c, its a linear equation
        if(a==0){
            return -c/b;
        }else{ //if a isnt zero, going to return the smallest t
            double delta=b*b-4*a*c;

            if(delta <0){ // if delta is smaller than 0, theres no solution, its imaginary
                // returns none, no solution
                return Double.NaN;
            }

            if(delta==0){ // if delta is zero only 1 solution
                return (-b)/(2*a);
            }

            t1=(-b+Math.sqrt(delta))/(2*a);
            t2=(-b-Math.sqrt(delta))/(2*a);

            //return the smallest of the time calculated but still bigger than zero
            if(t1 < 0){
                return t2;
            }

            if(t2 < 0){
                return t1;
            }

            //if both arent smaller than zero return the smallest
            return Math.min(t1, t2);
        }

    }


    // computing alpha
    public double computeAlpha(Particle A, Particle B, Particle C, double t){

        Vector2d Astar= new Vector2d();
        Vector2d Bstar= new Vector2d();
        Vector2d Cstar= new Vector2d();

        // adds t*velocity+p = new position of vector
        Astar.scaleAdd(t, A.v, A.p);
        Bstar.scaleAdd(t, B.v, B.p);
        Cstar.scaleAdd(t, C.v, C.p);

        Cstar.sub(Bstar); // Cstar now became C*-B*
        Astar.sub(Bstar); // Astar now became A*-B*

        //cstar dot became the top part of the equation, divided byt the lengthsquared of Astar-Bstar
        return Cstar.dot(Astar)/Astar.lengthSquared();


    }

    public Vector2d computeNormal(Particle A, Particle B, Particle C, double t){ // computes normal at the time of colision
        Vector2d Astar= new Vector2d();
        Vector2d Bstar= new Vector2d();

        // position at the time of colision
        Astar.scaleAdd(t, A.v, A.p);
        Bstar.scaleAdd(t, B.v, B.p);

        //subtraction of A and B, normalize it now
        Astar.sub(Bstar);

        Astar.normalize();

        //makes Bstar be the subtraction of position B and C
        // B - C
        Bstar.sub(B.p, C.p);

        //rotate it into 90 degress
        Astar.set(Astar.y*-1, Astar.x);

        // checks if it points towards C
        // if it doenst then flips it
        if(Astar.dot(Bstar) >0){
            Astar.scale(-1); // it flips the normal
        }

        return Astar;
    }

    public Double computeImpulse(Particle A, Particle B, Particle C, double alpha, Vector2d normal){ // computes impulse
        Vector2d relativeVelocity = computeRelVelocity(A, B, C, alpha);
        Double Vreln;

        //  normal.dot(relativeVelocity) = vreln-
        Vreln = normal.dot(relativeVelocity);

        // returns impulse
        return -(1+restitutionValue.getValue())*Vreln/(1/C.getMass()+alpha*alpha/A.getMass()+(1-alpha)*(1-alpha)/B.getMass());

    }

    public Vector2d computeRelVelocity(Particle A, Particle B, Particle C, double alpha){
        Vector2d vAlpha = new Vector2d();
        vAlpha.set(A.v);
        vAlpha.scale(alpha);
        // scale 1-alpha by B.v then add vAlpha then saves on vAlpha
        vAlpha.scaleAdd(1-alpha, B.v, vAlpha);
        vAlpha.sub(C.v, vAlpha);

        //relative velocity
        return vAlpha;
    }

    public void distributeImpulse(Particle A, Particle B, Particle C, double impulse, Vector2d directionImpulse, double alpha){
        // distributes the impulse for the new velocity, directionimpulse is the normal
        A.v.scaleAdd(-alpha*impulse/A.getMass(), directionImpulse, A.v);
        B.v.scaleAdd(-(1-alpha)*impulse/B.getMass(), directionImpulse, B.v);
        C.v.scaleAdd(impulse/C.getMass(), directionImpulse, C.v);
    }
    
}
