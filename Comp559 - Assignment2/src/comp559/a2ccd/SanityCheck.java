package comp559.a2ccd;

/**
 * Implementation of a forgiving edge-edge intersection test.  Guaranteed to return false
 * only if there are definitely two edges intersecting.
 * @author kry
 */
public class SanityCheck {

    private static double eps = 1e-4;
    
    /**
     * Performs an approximate sanity check.  This code errs on the side of forgiveness
     * so as not to stop the simulation inappropriately.
     * @param system 
     * @return true if no intersection, false if things "appear" to be OK
     */
    public static boolean sanityCheck( ParticleSystem system ) {
        for ( Spring s1 : system.springs ) {
            for ( Spring s2 : system.springs ) {
                if ( s1 == s2 ) continue;

                // discard the case where the two springs share a particle.
                // if they do happen to be superimposed, then there is likely to be
                // trouble at the next time step.  The one problem not caught will be
                // two segments passing through one another, but at least this leaves the 
                // simulation in a sane non-interpenetrating configuration.
                if ( s1.A == s2.A || s1.A == s2.B || s1.B == s2.A || s1.B == s2.B ) continue;

                // we'll solve a little linear system for the intersection of two lines, 
                // but give up if the determinant is too close to zero
                
                double a =  s1.B.p.x - s1.A.p.x;
                double b = -s2.B.p.x + s2.A.p.x;
                double c =  s1.B.p.y - s1.A.p.y;
                double d = -s2.B.p.y + s2.A.p.y;
                 
                double det = a*d-b*c;
                if ( Math.abs(det) < 1e-4 ) continue; // give up, but likely problems will be detected soon!
                double e = s2.A.p.x - s1.A.p.x;
                double f = s2.A.p.y - s1.A.p.y;
                double alpha1 = 1 / det * ( d*e - b*f);
                double alpha2 = 1 / det * (-c*e + a*f);
                
                if ( alpha1 > eps && alpha1 < 1-eps && alpha2 > eps && alpha2 < 1-eps ) {
                    // bad news!
                    return false;
                }
            }
        }
        return true;
    }
}
