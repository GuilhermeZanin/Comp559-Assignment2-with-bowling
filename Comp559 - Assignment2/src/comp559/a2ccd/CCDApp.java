package comp559.a2ccd;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.border.TitledBorder;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.util.gl2.GLUT;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.HorizontalFlowPanel;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;

/**
 * Template program for robust collision detection and response assignment
 * @author kry
 */
public class CCDApp implements SceneGraphNode, Interactor {
	
    private EasyViewer ev;
    
    public ParticleSystem system;

    public TestSystems testSystems;

    private double grabThresh = 10;

    public boolean sanity = true;
    
    /**
     * Entry point for application
     * @param args
     */
    public static void main(String[] args) {
        new CCDApp();        
    }
        
    /**
     * Creates the application / scene instance
     */
    public CCDApp() {
        system = new ParticleSystem();
        testSystems = new TestSystems( system );
        ev = new EasyViewer( "A2 - Robust Collision Processing", this, new Dimension(640,480), new Dimension(640,800) );
        ev.controlFrame.add("Test Systems and Factories", testSystems.getControls() );
        ev.controlFrame.add("Help", getHelpPanel() );
        ev.addInteractor(this);  // we add ourselves as an interactor to set up mouse and keyboard controls
    }
     
    @Override
    public void init(GLAutoDrawable drawable) {
        GL gl = drawable.getGL();
        gl.glEnable( GL.GL_BLEND );
        gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
        gl.glEnable( GL.GL_LINE_SMOOTH );
        gl.glEnable( GL2.GL_POINT_SMOOTH );
        system.init(drawable);
    }
        
    @Override
    public void display(GLAutoDrawable drawable) {

        GL2 gl = drawable.getGL().getGL2();
        if ( sanity ) {
        	gl.glClearColor( 1, 1, 1, 1);	
        } else {
        	gl.glClearColor( 1, 0.5f, 0.5f , 1);
        }
        gl.glClear( GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT );
        
        // First advance the system (if it is running or wants to be stepped)
        if ( run.getValue() || stepRequest ) {            
            // tell the testSystems about each step so that the alphabet soup factory
            // can generate new letters on the fly!
            testSystems.step();            
            for ( int i = 0; i < substeps.getValue(); i++ ) {
            	boolean resultOK = system.updateParticles( stepsize.getValue() / substeps.getValue() );
                if ( ! resultOK ) {
                    run.setValue( false );
                    sanity = false;
                    break;
                }
            }
            if ( ! SanityCheck.sanityCheck(system) ) {
                run.setValue( false );
                sanity = false;
            }            
        }
                    
        // We're doing 2D drawing only, so we'll call this helper function to set up a 2D projection 
        // where the origin is at the top left corner and the units are pixels
        EasyViewer.beginOverlay(drawable);
        
        gl.glTranslated(screenTranslationx.getValue(), screenTranslationy.getValue(), 0);
        
        gl.glDisable( GL2.GL_LIGHTING );
        
        system.display( drawable );

        // Here we'll display some extra stuff for the interface
        if ( system.grabbed ) { // draw the mouse spring
        	gl.glColor4d(0,0,1,0.75);
        	gl.glBegin( GL.GL_LINES);
            gl.glVertex2d( system.mouseSpring.A.p.x, system.mouseSpring.A.p.y );
            gl.glVertex2d( system.mouseSpring.B.p.x, system.mouseSpring.B.p.y );
            gl.glEnd();
        } else if ( mouseInWindow ) { // otherwise show the closest particle for grabbing
            findCloseParticles( xcurrent, ycurrent );
            if ( p1 != null && d1 < grabThresh ) {
                gl.glPointSize( 15f );
                gl.glColor4d(0,1,0,0.95);
                gl.glBegin( GL.GL_POINTS );
                gl.glVertex2d( p1.p.x, p1.p.y );
                gl.glEnd();        
            }
        }

        // Finally we'll display a string with useful information
        // about the system and the current stepping        
        if ( showInfo.getValue() ) {
        	// TODO: add your name to the display
            String text = "Guilherme Eidam Gaya Zanin\n" + system.toString();
        	text += "\n" + "h = " + stepsize.getValue() + "\n" +
                      "substeps = " + (int)(double) substeps.getValue() ;
	        gl.glColor3f(0.5f,0.5f,0.5f);
	        EasyViewer.printTextLines( drawable, text, 10, 15, 18, GLUT.BITMAP_9_BY_15 );
	        //EasyViewer.printTextLines( drawable, text );
        }

        EasyViewer.endOverlay(drawable);    
       
        // If we're recording, we'll save the step to an image file.
        // we'll also clear the step request here.
        if ( run.getValue() || stepRequest ) {
            stepRequest = false;        
            if ( record.getValue() ) {
                // write the frame
                File file = new File( "stills/" + dumpName + format.format(nextFrameNum) + ".png" );                                             
                nextFrameNum++;
                file = new File(file.getAbsolutePath().trim());
                ev.snapshot(drawable, file);
            }
        }
    }
    
    /** Flag to request the system be stepped in the next display callback */
    private boolean stepRequest = false;
        
    /** Base name of images to save */
    private String dumpName = "img";
    
    /** Index for the frame number we're saving. */
    private int nextFrameNum = 0;
    
    /** For formating the image file name when recording frames  */
    private NumberFormat format = new DecimalFormat("00000");
            
    private BooleanParameter record = new BooleanParameter( "record each step to image file", false );

    private BooleanParameter run = new BooleanParameter( "simulate (press SPACE in canvas to toggle)", false );
    
    private DoubleParameter stepsize = new DoubleParameter( "step size", 0.05, 1e-5, 1 );
    
    private IntParameter substeps = new IntParameter( "sub steps (integer)", 1, 1, 100);
        
    private BooleanParameter showInfo = new BooleanParameter( "display information overlay", true );
    
    private DoubleParameter screenTranslationx = new DoubleParameter( "screen translation x" , 0, -500,500 );
    
    private DoubleParameter screenTranslationy = new DoubleParameter( "screen translation y" , 0, -500,500 );
    
    public JPanel getHelpPanel() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        JTextArea ta = new JTextArea(
        		"   SPACE   - start/stop simulation \n" +
        		"   ENTER   - start/stop recording to stills \n" +
        		"   S       - step simulation \n" +
        		"   UP/DOWN - increase decrease sub-steps and step size \n" +
        		"   R       - reset to initial conditions \n" +
        		"   C       - clear all particles \n" +
        		"   Mouse Controls \n" +
        		"     - left drag for mouse spring on highlighted particle\n" +
        		"     - left click to pin/unpin \n" );
        ta.setEditable(false);
        vfp.add( ta );
        return vfp.getPanel();
    }
    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();

        HorizontalFlowPanel hfp = new HorizontalFlowPanel();
        JButton res1 = new JButton("640x360");
		hfp.add(res1);
		res1.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				ev.glCanvas.setSize(640, 360);
				ev.frame.setSize(ev.frame.getPreferredSize());
			}
		});
		JButton res2 = new JButton("1280x720");
		hfp.add(res2);
		res2.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				ev.glCanvas.setSize(1280, 720);
				ev.frame.setSize(ev.frame.getPreferredSize());
			}
		});
		hfp.add( record.getControls() );        
		vfp.add( hfp.getPanel() );
        
        VerticalFlowPanel vfp0 = new VerticalFlowPanel();
        vfp0.setBorder( new TitledBorder("Viewing Parameters" ) );
        ((TitledBorder) vfp0.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp0.add( showInfo.getControls() );
        vfp0.add( system.drawParticles.getControls() );
        vfp0.add( system.pointSize.getSliderControls(false) );
        vfp0.add( screenTranslationx.getSliderControls(false) );
        vfp0.add( screenTranslationy.getSliderControls(false) );
        vfp.add( vfp0.getPanel() );
        
        VerticalFlowPanel vfp1 = new VerticalFlowPanel();
        vfp1.setBorder( new TitledBorder("Numerical Integration Controls"));
        ((TitledBorder) vfp1.getPanel().getBorder()).setTitleFont(new Font("Tahoma", Font.BOLD, 18));
        vfp1.add( run.getControls() );        
        vfp1.add( stepsize.getSliderControls(true) );
        vfp1.add( substeps.getControls() );
        vfp.add( vfp1.getPanel() );
        
        vfp.add( system.getControls() );
               
        return vfp.getPanel();
    }
    
    // Some member variables to help us keep track of close particles
    private Particle p1 = null;
    private Particle p2 = null;
    private double d1 = 0;
    private double d2 = 0;
    
    /**
     * Finds the two closest particles for showing potential spring connections
     * @param x 
     * @param y 
     */
    private void findCloseParticles( int x, int y ) {
        List<Particle> particles = system.particles;
        p1 = null;
        p2 = null;
        d1 = 0;
        d2 = 0;
        if ( particles.size() > 0 ) {
            for ( Particle p : particles ) {                
                double d = p.distance( x, y );
                if ( p1 == null || d < d1 ) {
                    p2 = p1; d2 = d1; p1 = p; d1 = d;
                } else if ( p2 == null || d < d2 ) {
                    p2 = p; d2 = d;
                }
            }      
        }
    }  
    
    private int xcurrent = 0;
    private int ycurrent = 0;
    private boolean mouseInWindow = false;
    
    @Override
    public void attach(Component component) {
        component.addMouseMotionListener( new MouseMotionListener() {
            @Override
            public void mouseDragged(MouseEvent e) {
            	system.mouseParticle.p.set( e.getPoint().x, e.getPoint().y );
        		xcurrent = e.getPoint().x;
        		ycurrent = e.getPoint().y;
            }
            @Override
            public void mouseMoved(MouseEvent e) {
                xcurrent = e.getPoint().x;
                ycurrent = e.getPoint().y;
            }
        } );
        component.addMouseListener( new MouseListener() {
            @Override
            public void mouseClicked(MouseEvent e) {
            	if ( e.getButton() == MouseEvent.BUTTON1 ) {
            		findCloseParticles( e.getPoint().x, e.getPoint().y );
                    if ( p1 != null && d1 < grabThresh ) {
                        p1.pinned = !p1.pinned;
            		}
            	}
            }
            @Override
            public void mouseEntered(MouseEvent e) {
                mouseInWindow = true;
            }
            @Override
            public void mouseExited(MouseEvent e) {
                mouseInWindow = false;
            }
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.getButton() == MouseEvent.BUTTON1 ) {
	                int x = e.getPoint().x;
	                int y = e.getPoint().y;
	                findCloseParticles( x, y );
	                if ( p1 != null && d1 < grabThresh ) {
	                	system.mouseSpring.B = p1;
	                	system.mouseParticle.p.set( x, y );       
	                	system.grabbed = true;
	                }
                }
            }
            @Override
            public void mouseReleased(MouseEvent e) {
                if (e.getButton() == MouseEvent.BUTTON1 ) {
                	system.mouseSpring.B = null;
                	system.grabbed = false;
                }
            }
        } );
        component.addKeyListener( new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
                    run.setValue( ! run.getValue() ); 
                } else if ( e.getKeyCode() == KeyEvent.VK_S ) {                    
                    stepRequest = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_R ) {
                    system.resetParticles();
                    sanity = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_C ) {                   
                    system.clear();
                    testSystems.resetFactory();
                    p1 = null;
                    p2 = null;
                    sanity = true;
                } else if ( e.getKeyCode() == KeyEvent.VK_UP ) {
                    int ss = substeps.getValue();
                    if ( ss == substeps.getMaximum() ) return;
                    substeps.setValue( ss + 1 );
                    stepsize.setValue( stepsize.getValue() * (ss+1)/ss );
                } else if ( e.getKeyCode() == KeyEvent.VK_DOWN ) {
                    int ss = substeps.getValue();
                    if ( ss == substeps.getMinimum() ) return;
                    substeps.setValue( ss - 1 );
                    stepsize.setValue( stepsize.getValue() *(ss-1)/ss );
                } else if ( e.getKeyCode() == KeyEvent.VK_ESCAPE ) {
                    // quit the program
                    ev.stop();
                } else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
                    // toggle recording of steps to png files
                    record.setValue( ! record.getValue() );
                }
            }
        } );
    }
    
}
