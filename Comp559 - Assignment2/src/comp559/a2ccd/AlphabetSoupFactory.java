package comp559.a2ccd;

import java.awt.Font;
import java.awt.GraphicsEnvironment;
import java.awt.Shape;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.font.FontRenderContext;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.PathIterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.border.TitledBorder;
import javax.vecmath.Point2d;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.parameters.Parameter;
import mintools.parameters.ParameterListener;
import mintools.swing.VerticalFlowPanel;

/**
 * Letter generation helper class.  Deals with ugly issues such as trying to get an 
 * approximately constant segment length piece-wise representation of letters on both 
 * straight and curved segements.  It also deals with closing loops by searching for
 * previously generated particles that are within a threshold of a particle to be created.
 * 
 * Note that the constant size segment length is likely to make for very challenging collision
 * response problems!
 *  
 * @author kry
 */
public class AlphabetSoupFactory {

    IntParameter which = new IntParameter( "which font", 12, 0, 73);
    DoubleParameter size = new DoubleParameter( "size", 200, 10, 1000 );
    //DoubleParameter curveResolution = new DoubleParameter( "curve resolution", 1, 1, 20);
    BooleanParameter bold = new BooleanParameter( "bold", false );
    BooleanParameter italic = new BooleanParameter( "italic", false );
    
    DoubleParameter offsetx = new DoubleParameter("x offset", 50, 0, 400 );
    DoubleParameter offsety = new DoubleParameter("y offset", 250, 0, 400 );
    
    IntParameter segmentBreak = new IntParameter("approximate segment length", 10, 5, 100 );
    
    JTextField text = new JTextField("test");
    JLabel fname = new JLabel();
    
    String[] fontNames;
    
    ParticleSystem system;
    
    /**
     * Creates a new alphabet soup factory
     */
    public AlphabetSoupFactory( ParticleSystem system ) {
        // Get the local graphics environment
        GraphicsEnvironment env = GraphicsEnvironment.getLocalGraphicsEnvironment();      
        //Get the font names from the graphics environment
        fontNames = env.getAvailableFontFamilyNames();        
        which.setMaximum(fontNames.length -1);
        fname.setText(fontNames[(int)which.getValue()]);
        which.addParameterListener(new ParameterListener<Integer>() {
           @Override
           public void parameterChanged(Parameter<Integer> parameter) {
               fname.setText(fontNames[(int)which.getValue()]);   
           }        
        });
        this.system = system;
    }
    
    /**
     * Gets the letter creation controls
     * @return the control panel
     */
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.setBorder( new TitledBorder("Letter Generation Parameters") );
        ((TitledBorder) vfp.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp.add( fname );
        vfp.add( which.getSliderControls() );        
        vfp.add( size.getSliderControls(true));
        //vfp.add( curveResolution.getSliderControls(false));
        vfp.add( bold.getControls() );
        vfp.add( italic.getControls() );
        vfp.add( text );
        
        vfp.add( segmentBreak.getSliderControls() );
        vfp.add( offsetx.getSliderControls(false) );
        vfp.add( offsety.getSliderControls(false) );
        
        JButton createLetters = new JButton( "Create Letters" );
        createLetters.addActionListener( new ActionListener() {			
			@Override
			public void actionPerformed(ActionEvent e) {
				createLetter(system);
		        system.name = "test letters";				
			}
		});
        vfp.add(createLetters);
        return vfp.getPanel();   
    }
    
    /** 
     * Creates the letters specified in the control panel
     * @param system
     */
    public void createLetter( ParticleSystem system ) {
        createLetter( system, text.getText(), offsetx.getValue(), offsety.getValue() );
    }
    
    /**
     * Creates the letter specified in the control panel
     * @param system
     * @param letters 
     * @param offx 
     * @param offy 
     */
    public void createLetter( ParticleSystem system, String letters, double offx, double offy ) {
                
        // perhaps not too efficient...
        
        // we'll get the stroke information for the given text and the given 
        // font on every display call
        int style = 0;
        style |= bold.getValue() ? Font.BOLD : 0;
        style |= italic.getValue() ? Font.ITALIC : 0;                       
        fname.setText(fontNames[(int)which.getValue()]);        
        Font f = new Font( fname.getText(), style, (int)(double)size.getValue() );        
        FontRenderContext defaultFRC = new FontRenderContext(null, false, false);
        GlyphVector gv = f.createGlyphVector(defaultFRC, letters );        
        Shape shape = gv.getOutline();
        AffineTransform at = new AffineTransform();
        at.setToIdentity();
        PathIterator pit = shape.getPathIterator( at );
        
        p2 = null;
        p1 = null;
        p0 = null;
        newParticles.clear();
        
        int segmentBreakSize = segmentBreak.getValue();
        
        // int n = (int) curveResolution.getValue();
        double[] c = new double[6];
        Point2d p = new Point2d(0,0);                
        while ( !pit.isDone() ) {            
            int type = pit.currentSegment(c);
            if ( type == PathIterator.SEG_LINETO ) {
                
                
                Point2d nextp = new Point2d( c[0], c[1] );
                double dist = p.distance(nextp);
                if ( dist > segmentBreakSize ) {
                    int N = (int)dist / segmentBreakSize;
                    for (int i = 1; i < N; i++ ) {
                        Point2d tmp = new Point2d();
                        double alpha = (double)i/(double)N;
                        tmp.interpolate( nextp, p, 1- alpha);
                        addParticle( tmp.x, tmp.y, offx, offy, system );
                    }
                }
                
                addParticle( c[0], c[1], offx, offy, system );
                p.set( c[0], c[1]);                
            } else if ( type == PathIterator.SEG_MOVETO ) {
                p0 = null;
                p1 = null;
                p2 = null;
                addParticle( c[0], c[1], offx, offy, system );
                p.set( c[0], c[1]);                
            } else if ( type == PathIterator.SEG_CLOSE ) {
                p0 = null;
                p1 = null;
                p2 = null;
                p.set( c[0], c[1]);
            } else if ( type == PathIterator.SEG_CUBICTO ) {
                //  P(t) = B(3,0)*CP + B(3,1)*P1 + B(3,2)*P2 + B(3,3)*P3,  0 <= t <= 1
                //         B(n,m) = mth coefficient of nth degree Bernstein polynomial
                //                = C(n,m) * t^(m) * (1 - t)^(n-m)
                //         C(n,m) = Combinations of n things, taken m at a time
                //                = n! / (m! * (n-m)!)
                
                double x = 0;
                double y = 0;
                
                double dist = approximateCubicCurveLength(p, c);
                int N = (int) (dist / segmentBreakSize);
                if ( N < 1 ) N = 1;
                
                for ( int i = 1; i <= N; i++) { // don't start at the beginning... ?  i.e., i = 1, not zero
                    double t = i/(double)N;                    
                    double b0 = 1 * t*t*t;
                    double b1 = 3 * t*t*(1-t);
                    double b2 = 3 * t*(1-t)*(1-t);
                    double b3 = 1 * (1-t)*(1-t)*(1-t);                    
                    x = p.x*b3 + c[0]*b2 + c[2]*b1 + c[4]*b0;
                    y = p.y*b3 + c[1]*b2 + c[3]*b1 + c[5]*b0;
                    addParticle( x, y, offx, offy, system );
                }
                p.set( x, y );                
            } else if ( type == PathIterator.SEG_QUADTO ) {
                // P(t) = B(2,0)*CP + B(2,1)*P1 + B(2,2)*P2, 0 <= t <= 1
                //        B(n,m) = mth coefficient of nth degree Bernstein polynomial
                //               = C(n,m) * t^(m) * (1 - t)^(n-m)
                //        C(n,m) = Combinations of n things, taken m at a time
                //               = n! / (m! * (n-m)!)

                double dist = approximateQuadraticCurveLength(p, c);
                int N = (int) (dist / segmentBreakSize);
                if ( N < 1 ) N = 1;
                
                double x = 0;
                double y = 0;
                for ( int i = 1; i <= N; i++) { // don't start at the beginning!
                    double t = i/(double)N;                    
                    double b0 = 1 * t*t;
                    double b1 = 2 * t*(1-t);
                    double b2 = 1 * (1-t)*(1-t);                                       
                    x = p.x*b2 + c[0]*b1 + c[2]*b0;
                    y = p.y*b2 + c[1]*b1 + c[3]*b0;
                    addParticle( x, y, offx, offy, system );
                }
                p.set( x, y );                
            }
            pit.next();
        }
    }   
    
    private double approximateCubicCurveLength( Point2d p, double[] c ) {
        Point2d tmp1 = new Point2d();
        Point2d tmp2 = new Point2d();
        double length = 0;
        // 5 segments should be plenty for approximating the length.
        tmp1.set( p );
        for ( int i = 1; i <= 5; i++) { 
            double t = i/(double)5;                    
            double b0 = 1 * t*t*t;
            double b1 = 3 * t*t*(1-t);
            double b2 = 3 * t*(1-t)*(1-t);
            double b3 = 1 * (1-t)*(1-t)*(1-t);                    
            double x = p.x*b3 + c[0]*b2 + c[2]*b1 + c[4]*b0;
            double y = p.y*b3 + c[1]*b2 + c[3]*b1 + c[5]*b0;
            tmp2.set( x, y );            
            length += tmp1.distance(tmp2);
            tmp1.set( tmp2 );
        }
        return length;
    }
    
    private double approximateQuadraticCurveLength( Point2d p, double[] c ) {
        Point2d tmp1 = new Point2d();
        Point2d tmp2 = new Point2d();
        double length = 0;
        // 5 segments should be plenty for approximating the length.
        tmp1.set( p );
        for ( int i = 1; i <= 5; i++) { 
            double t = i/(double)5;                    
            double b0 = 1 * t*t;
            double b1 = 2 * t*(1-t);
            double b2 = 1 * (1-t)*(1-t);                                       
            double x = p.x*b2 + c[0]*b1 + c[2]*b0;
            double y = p.y*b2 + c[1]*b1 + c[3]*b0;
            tmp2.set( x, y );            
            length += tmp1.distance(tmp2);
            tmp1.set( tmp2 );
        }
        return length;
    }
    
    Particle p2, p1, p0;

    List<Particle> newParticles = new LinkedList<Particle>();
    
    /**
     * Creates a new particle if necessary and joins it to the previous with springs and leaf springs
     * @param x
     * @param y
     * @param system
     */
    private void addParticle( double x, double y, double offx, double offy, ParticleSystem system ) {        
        // gross! sloth! but needed a solution quickly... 
        Particle newp2 = findCloseParticle( x + offx, y + offy );
        if ( newp2 == p2 && p2 != null ) {
            return; // skip it... it was too close to the last!            
        }              
        p0 = p1;
        p1 = p2;
        p2 = newp2;
        if ( p2 == null ) {
            p2 = new Particle( x + offx, y + offy, 0, 0 );
            newParticles.add( p2 );
            system.particles.add( p2 );
        }
        if ( p1 != null ) system.springs.add( new Spring( p1, p2 ) );
        if ( p0 == p1 && p0 != null) {
            return;
        }
        if ( p0 != null ) system.bendingSprings.add( new BendingSpring( p0, p1, p2 ) );
    }
    
    private Particle findCloseParticle( double x, double y ) {
        Point2d a = new Point2d(x,y);
        double minDist = Double.POSITIVE_INFINITY;
        Particle bestParticle = null;
        for ( Particle p : newParticles ) {
            double d = a.distance( p.p );
            if ( d < 4 && d < minDist ) {
                minDist = d;
                bestParticle = p;
            }            
        }
        return bestParticle;
    }
    
}
