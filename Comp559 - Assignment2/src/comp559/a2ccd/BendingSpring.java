package comp559.a2ccd;

/**
 * Leaf Spring class for COMP 599 assignment 2
 * (do not distribute)
 * @author kry
 */
public class BendingSpring {

    Particle A = null;
    Particle B = null;
    Particle C = null;    
    
    /** Modulation for bending stiffness */
    static double kbMod = 1e4;

    /** Modulation for bending damping */
    static double kbdMod = 1;
    
    /** Spring bending base stiffness */
    double kbBase = 1;
    
    /** Spring bending stiffness */
    double kb = kbBase * kbMod;
    
    /** Spring damping coefficient */
    double kbd = 1;
    
    /** Rest angle */
    double theta0 = 0;
    
    /** Rest curvature */
    double kappa0 = 0;
        
    /**
     * Creates a spring connecting two particles.
     * The rest length should be set
     * @param p0
     * @param p1
     * @param p2
     */
    public BendingSpring( Particle p0, Particle p1, Particle p2 ) {
        this.A = p0;
        this.B = p1;
        this.C = p2;
        computeRestAngle();
        computeRestCurvature();
    }
    
    /**
     * Compute the rest curvature of the connected particles
     */
    public void computeRestCurvature() {
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t9 = Bx * Bx;
    	double t12 = Ax * Ax;
    	double t13 = By * By;
    	double t16 = Ay * Ay;
    	double t18 = Math.sqrt(t9 - 2 * Bx * Ax + t12 + t13 - 2 * By * Ay + t16);
    	double t21 = Cx * Cx;
    	double t24 = Cy * Cy;
    	double t28 = Math.sqrt(t21 - 2 * Cx * Bx + t9 + t24 - 2 * Cy * By + t13);
    	double t35 = Math.sqrt(t21 - 2 * Cx * Ax + t12 + t24 - 2 * Cy * Ay + t16);
    	kappa0 = (2 * (Bx - Ax) * (Ay - Cy) + 2 * (By - Ay) * (Cx - Ax)) / t18 / t28 / t35;
    }
    
    /**
     * Computes the rest length of the connected particles
     */
    public void computeRestAngle() {
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t1 = Ay - By;
    	double t2 = Cx - Bx;
    	double t4 = Bx - Ax;
    	double t5 = Cy - By;
    	theta0 = Math.atan2(t1 * t2 + t4 * t5, t4 * t2 - t1 * t5);
    }
    
    public void apply( boolean useCurvature ) {
    	if ( useCurvature ) {
    		applyCurvatureBasedForce();
    	} else {
    		applyAngleBasedForce();
    	}
    }
    
    public void applyCurvatureBasedForce() {
    	kb = kbBase * kbMod;
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t1 = Bx - Ax;
    	double t2 = Ay - Cy;
    	double t3 = By - Ay;
    	double t4 = Cx - Ax;
    	double t5 = 2 * t1 * t2 + 2 * t3 * t4;
    	double t6 = Bx * Bx;
    	double t7 = Ax * Ax;
    	double t8 = By * By;
    	double t9 = Ay * Ay;
    	double t10 = 2 * By;
    	double t11 = 2 * Bx;
    	double t12 = t6 - t11 * Ax + t7 + t8 - t10 * Ay + t9;
    	double t13 = Math.pow(t12, -0.3e1 / 0.2e1);
    	 t12 = t12 * t13;
    	double t14 = Cy * Cy;
    	double t15 = Cx * Cx;
    	t6 = t15 - t11 * Cx + t6 + t14 - t10 * Cy + t8;
    	t8 = Math.pow(t6, -0.3e1 / 0.2e1);
    	t6 = t6 * t8;
    	t7 = t15 - 2 * Cx * Ax + t7 + t14 - 2 * Cy * Ay + t9;
    	t9 = Math.pow(t7, -0.3e1 / 0.2e1);
    	t7 = t7 * t9;
    	t10 = 2 * Cy - 2 * By;
    	t1 = 2 * t1;
    	t11 = -t1;
    	t4 = 2 * t4;
    	t14 = t12 * t9;
    	t15 = t13 * t7;
    	double t16 = -2 * Cx + 2 * Bx;
    	t3 = 2 * t3;
    	t2 = 2 * t2;
    	double t17 = t16 * t12;
    	t13 = t13 * t6;
    	t9 = t6 * t9;
    	double t18 = t8 * t7;
    	double t19 = kb * (t5 * t12 * t6 * t7 - kappa0);
    	A.f.x += -t19 * t6 * ((-t15 * t11 / 2 + t14 * t4 / 2) * t5 + t10 * t12 * t7);
    	A.f.y += -t19 * t6 * ((t15 * t3 / 2 - t14 * t2 / 2) * t5 + t17 * t7);
    	B.f.x += -t19 * t7 * ((-t13 * t1 / 2 - t17 * t8 / 2) * t5 + t2 * t12 * t6);
    	B.f.y += -t19 * t7 * ((-t13 * t3 / 2 + t12 * t8 * t10 / 2) * t5 + t4 * t12 * t6);
    	C.f.x += -t19 * t12 * ((t18 * t16 / 2 - t9 * t4 / 2) * t5 + t3 * t6 * t7);
    	C.f.y += -t19 * t12 * ((-t18 * t10 / 2 + t9 * t2 / 2) * t5 + t11 * t6 * t7);

    }
    
    public void applyAngleBasedForce() {
    	kb = kbBase * kbMod;
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t1 = Ay - By;
    	double t2 = Cx - Bx;
    	double t3 = Bx - Ax;
    	double t4 = Cy - By;
    	double t5 = t1 * t2 + t3 * t4;
    	double t6 = -t1;
    	double t7 = t3 * t2 + t6 * t4;
    	t4 = -t4;
    	double t8 = 1 / t7;
    	double t9 = t8 * t5;
    	double t10 = 1 + t5*t5 * t8*t8;
    	t10 = 0.1e1 / t10;
    	t5 = kb * (Math.atan2(t5, t7) - theta0);
    	if ( Math.abs(t5) < 1e-5 ) return;
    	// note t8 is probably infinity in this case?
    	A.f.x += -t5 * t8 * (t4 + t9 * t2) * t10;
    	A.f.y += -t5 * t8 * (t2 - t9 * t4) * t10;
    	B.f.x += -t5 * t8 * (Cy - Ay - t9 * (Cx - 2 * Bx + Ax)) * t10;
    	B.f.y += -t5 * t8 * (-Cx + Ax - t9 * (Cy - 2 * By + Ay)) * t10;
    	C.f.x += -t5 * t8 * (t1 - t9 * t3) * t10;
    	C.f.y += -t5 * t8 * (t3 - t9 * t6) * t10;
    }
    
    /** Stiffness matrix */
    private double[][] K = new double[6][6];
    
    public void computeAngleBasedStiffnessKdx() {
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t1 = -Cy + By;
    	double t2 = Bx - Ax;
    	double t3 = Cx - Bx;
    	double t4 = By - Ay;
    	double t5 = -t1;
    	double t6 = t2 * t3 + t4 * t5;
    	t5 = -t3 * t4 + t2 * t5;
    	double t7 = -t3;
    	double t8 = 0.1e1 / t6;
    	double t9 = t5 * t8;
    	double t10 = t9 * t7;
    	double t11 = t1 - t10;
    	double t12 = t8 * t11;
    	double t13 = Math.pow(t8, 2);
    	double t14 = 1 + Math.pow(t5, 2) * t13;
    	t6 = Math.atan2(t5, t6) - theta0;
    	double t15 = 2 * t13;
    	double t16 = t15 * t5;
    	t14 = 0.1e1 / t14;
    	double t17 = t14 * kb;
    	double t18 = t9 * t1;
    	double t19 = t3 - t18;
    	double t20 = t8 * t19;
    	t18 = -t3 + 2 * t18;
    	double t21 = t6 * t16 * t19;
    	double t22 = t17 * ((-t20 + t21) * t14 * t12 - t6 * t13 * (t18 * t7 - Math.pow(t1, 2)));
    	double t23 = Cy - Ay;
    	double t24 = Cx - 2 * Bx + Ax;
    	double t25 = t9 * t24;
    	double t26 = t23 - t25;
    	double t27 = t8 * t26;
    	t10 = -t1 + 2 * t10;
    	double t28 = t6 * t16 * t26;
    	double t29 = (-t27 + t28) * t14;
    	double t30 = t17 * (t29 * t12 - t6 * t13 * (t10 * t24 - t23 * t7 - t5));
    	double t31 = -Cx + Ax;
    	double t32 = Cy - 2 * By + Ay;
    	double t33 = t31 - t9 * t32;
    	double t34 = t8 * t33;
    	double t35 = t6 * t16 * t33;
    	double t36 = (-t34 + t35) * t14;
    	t10 = t17 * (t36 * t12 - t6 * t8 * (t8 * (t10 * t32 - t31 * t7) + 1));
    	double t37 = t9 * t2;
    	double t38 = -t4 - t37;
    	double t39 = t8 * t38;
    	t37 = t4 + 2 * t37;
    	double t40 = t1 * t2;
    	double t41 = t6 * t16 * t38;
    	double t42 = (-t39 + t41) * t14;
    	double t43 = t17 * (t42 * t12 - t6 * t13 * (t37 * t7 - t40 + t5));
    	t9 = t9 * t4;
    	double t44 = t2 - t9;
    	double t45 = t8 * t44;
    	t9 = -t2 + 2 * t9;
    	double t46 = t1 * t4;
    	double t47 = t6 * t16 * t44;
    	double t48 = (-t45 + t47) * t14;
    	double t49 = t17 * (t48 * t12 - t6 * t8 * (t8 * (t9 * t7 - t46) - 1));
    	t29 = t17 * (t29 * t20 - t6 * t8 * (t8 * (t18 * t24 - t23 * t1) - 1));
    	double t50 = t17 * (t36 * t20 - t6 * t13 * (t18 * t32 - t31 * t1 - t5));
    	t18 = t17 * (t42 * t20 - t6 * t8 * (t8 * (t18 * t2 + t1 * t4) + 1));
    	t3 = t17 * (t48 * t20 - t6 * t13 * ((2 * t46 * t8 + 1) * t5 - t3 * t4 - t40));
    	t25 = t17 * (t36 * t27 - t6 * t13 * ((-t23 + 2 * t25) * t32 - t31 * t24));
    	t36 = t17 * (t42 * t27 - t6 * t13 * (t37 * t24 - t2 * t23 - t5));
    	t23 = t17 * (t48 * t27 - t6 * t8 * (t8 * (t9 * t24 - t23 * t4) + 1));
    	t8 = t17 * (t42 * t34 - t6 * t8 * (t8 * (t37 * t32 - t2 * t31) - 1));
    	t9 = t17 * (t48 * t34 - t6 * t13 * (t9 * t32 - t31 * t4 - t5));
    	t13 = t17 * (t48 * t39 - t6 * t13 * (t37 * t4 - t2*t2 ));
    	K[0][0] = t17 * (-Math.pow(t12, 2) * t14 + t6 * t15 * t7 * t11 + t6 * t12 * t14 * t16 * t11);
    	K[0][1] = t22;
    	K[0][2] = t30;
    	K[0][3] = t10;
    	K[0][4] = t43;
    	K[0][5] = t49;
    	K[1][0] = t22;
    	K[1][1] = t17 * (-Math.pow(t20, 2) * t14 + t6 * t15 * t1 * t19 + t21 * t20 * t14);
    	K[1][2] = t29;
    	K[1][3] = t50;
    	K[1][4] = t18;
    	K[1][5] = t3;
    	K[2][0] = t30;
    	K[2][1] = t29;
    	K[2][2] = t17 * (-Math.pow(t27, 2) * t14 - t6 * t15 * (t5 - t26 * t24) + t28 * t27 * t14);
    	K[2][3] = t25;
    	K[2][4] = t36;
    	K[2][5] = t23;
    	K[3][0] = t10;
    	K[3][1] = t50;
    	K[3][2] = t25;
    	K[3][3] = t17 * (-Math.pow(t34, 2) * t14 - t6 * t15 * (t5 - t33 * t32) + t35 * t34 * t14);
    	K[3][4] = t8;
    	K[3][5] = t9;
    	K[4][0] = t43;
    	K[4][1] = t18;
    	K[4][2] = t36;
    	K[4][3] = t8;
    	K[4][4] = t17 * (-Math.pow(t39, 2) * t14 + t6 * t15 * t2 * t38 + t41 * t39 * t14);
    	K[4][5] = t13;
    	K[5][0] = t49;
    	K[5][1] = t3;
    	K[5][2] = t23;
    	K[5][3] = t9;
    	K[5][4] = t13;
    	K[5][5] = t17 * (-Math.pow(t45, 2) * t14 + t6 * t15 * t4 * t44 + t47 * t45 * t14);
    }
    
    public void computeCurvatureBasedStiffnessKdx() {
    	double Ax = A.p.x;
    	double Ay = A.p.y;
    	double Bx = B.p.x;
    	double By = B.p.y;
    	double Cx = C.p.x;
    	double Cy = C.p.y;
    	double t1 = Cy - By;
    	double t2 = Bx * Bx;
    	double t5 = Ax * Ax;
    	double t6 = By * By;
    	double t9 = Ay * Ay;
    	double t10 = t2 - 2 * Bx * Ax + t5 + t6 - 2 * By * Ay + t9;
    	double t11 = Math.sqrt(t10);
    	double t12 = 0.1e1 / t11;
    	double t13 = 2 * t12 * t1;
    	double t14 = Cx * Cx;
    	double t17 = Cy * Cy;
    	double t20 = t14 - 2 * Cx * Bx + t2 + t17 - 2 * Cy * By + t6;
    	double t21 = Math.sqrt(t20);
    	double t22 = 0.1e1 / t21;
    	double t27 = t14 - 2 * Cx * Ax + t5 + t17 - 2 * Cy * Ay + t9;
    	double t28 = Math.sqrt(t27);
    	double t29 = 0.1e1 / t28;
    	double t30 = t29 * t22;
    	double t32 = Bx - Ax;
    	double t33 = Ay - Cy;
    	double t35 = By - Ay;
    	double t36 = Cx - Ax;
    	double t38 = t33 * t32 + t36 * t35;
    	double t40 = 0.1e1 / t11 / t10;
    	double t41 = 2 * t40 * t38;
    	double t42 = -2 * t32 * t30;
    	double t45 = 2 * t12 * t38;
    	double t47 = 0.1e1 / t28 / t27;
    	double t48 = t47 * t22;
    	double t49 = -2 * t36 * t48;
    	double t52 = t30 * t13 - t42 * t41 / 2 - t49 * t45 / 2;
    	double t53 = t52 * t52;
    	double t57 = (t30 * t45 - kappa0) * kb;
    	double t58 = 2 * t40 * t1;
    	double t61 = t10 * t10;
    	double t64 = 2 / t11 / t61 * t38;
    	double t65 = 4 * t32 * t32;
    	double t68 = 0.3e1 / 0.4e1 * t65 * t30 * t64;
    	double t69 = t22 * t41;
    	double t70 = -2 * t32 * t47;
    	double t74 = t30 * t41;
    	double t75 = t27 * t27;
    	double t77 = 0.1e1 / t28 / t75;
    	double t78 = t77 * t22;
    	double t79 = 4 * t36 * t36;
    	double t82 = 0.3e1 / 0.4e1 * t79 * t78 * t45;
    	double t83 = t48 * t45;
    	double t87 = -Cx + Bx;
    	double t88 = 2 * t12 * t87;
    	double t90 = -2 * t35 * t30;
    	double t93 = 2 * t33 * t48;
    	double t96 = t30 * t88 - t90 * t41 / 2 - t93 * t45 / 2;
    	double t103 = 2 * t40 * t87;
    	double t106 = t22 * t64;
    	double t107 = -2 * t32 * t29;
    	double t116 = -2 * t36 * t47;
    	double t120 = t22 * t45;
    	double t121 = -2 * t36 * t77;
    	double t127 = -t52 * t96 * kb - (-t90 * t58 / 2 - t93 * t13 / 2 - t42 * t103 / 2 - 0.3e1 / 0.2e1 * t35 * t107 * t106 + t33 * t70 * t69 / 2 - t49 * t88 / 2 - t35 * t116 * t69 / 2 + 0.3e1 / 0.2e1 * t33 * t121 * t120) * t57;
    	double t128 = 2 * t12 * t33;
    	double t130 = 2 * t32 * t30;
    	double t134 = 0.1e1 / t21 / t20;
    	double t135 = t29 * t134;
    	double t136 = 2 * t87 * t135;
    	double t139 = t30 * t128 - t130 * t41 / 2 - t136 * t45 / 2;
    	double t140 = t139 * kb;
    	double t145 = t136 * t13 / 2;
    	double t146 = 2 * t40 * t33;
    	double t152 = t134 * t41;
    	double t161 = t134 * t45;
    	double t167 = -t52 * t140 - (-t130 * t58 / 2 - t145 - t42 * t146 / 2 + 0.3e1 / 0.2e1 * t32 * t107 * t106 + t87 * t107 * t152 / 2 + t74 - t49 * t128 / 2 + t32 * t116 * t69 / 2 + t87 * t116 * t161 / 2) * t57;
    	double t168 = 2 * t12 * t36;
    	double t170 = 2 * t35 * t30;
    	double t173 = -2 * t1 * t135;
    	double t176 = t30 * t168 - t170 * t41 / 2 - t173 * t45 / 2;
    	double t177 = t176 * kb;
    	double t181 = 2 * t29 * t22 * t12;
    	double t186 = 2 * t40 * t36;
    	double t205 = -t52 * t177 - (-t181 - t170 * t58 / 2 - t173 * t13 / 2 - t42 * t186 / 2 + 0.3e1 / 0.2e1 * t35 * t107 * t106 - t1 * t107 * t152 / 2 - t49 * t168 / 2 + t35 * t116 * t69 / 2 - t1 * t116 * t161 / 2) * t57;
    	double t206 = 2 * t12 * t35;
    	double t208 = -2 * t87 * t135;
    	double t211 = 2 * t36 * t48;
    	double t214 = t30 * t206 - t208 * t45 / 2 - t211 * t45 / 2;
    	double t215 = t214 * kb;
    	double t221 = 2 * t40 * t35;
    	double t223 = t42 * t221 / 2;
    	double t240 = -t52 * t215 - (-t208 * t13 / 2 - t211 * t13 / 2 - t223 - t87 * t107 * t152 / 2 + t36 * t70 * t69 / 2 - t49 * t206 / 2 - t87 * t116 * t161 / 2 + 0.3e1 / 0.2e1 * t36 * t121 * t120 + t83) * t57;
    	double t241 = -2 * t12 * t32;
    	double t243 = 2 * t1 * t135;
    	double t246 = -2 * t33 * t48;
    	double t249 = t30 * t241 - t243 * t45 / 2 - t246 * t45 / 2;
    	double t250 = t249 * kb;
    	double t252 = 4 * t1 * t1;
    	double t277 = -t52 * t250 - (t181 - t135 * t12 * t252 / 2 - t246 * t13 / 2 - t30 * t40 * t65 / 2 + t1 * t107 * t152 / 2 - t33 * t70 * t69 / 2 - t49 * t241 / 2 + t1 * t116 * t161 / 2 - 0.3e1 / 0.2e1 * t33 * t121 * t120) * t57;
    	double t278 = t96 * t96;
    	double t282 = 4 * t35 * t35;
    	double t285 = 0.3e1 / 0.4e1 * t282 * t30 * t64;
    	double t286 = -2 * t35 * t47;
    	double t290 = 4 * t33 * t33;
    	double t293 = 0.3e1 / 0.4e1 * t290 * t78 * t45;
    	double t300 = 4 * t87 * t87;
    	double t306 = -2 * t35 * t29;
    	double t316 = 2 * t33 * t47;
    	double t325 = -t96 * t140 - (t181 - t130 * t103 / 2 - t135 * t12 * t300 / 2 - t90 * t146 / 2 + 0.3e1 / 0.2e1 * t32 * t306 * t106 + t87 * t306 * t152 / 2 - t48 * t12 * t290 / 2 + t32 * t316 * t69 / 2 + t87 * t316 * t161 / 2) * t57;
    	double t340 = t93 * t168 / 2;
    	double t349 = -t96 * t177 - (-t170 * t103 / 2 - t173 * t88 / 2 - t90 * t186 / 2 + 0.3e1 / 0.2e1 * t35 * t306 * t106 - t1 * t306 * t152 / 2 + t74 - t340 + t35 * t316 * t69 / 2 - t1 * t316 * t161 / 2) * t57;
    	double t368 = 2 * t33 * t77;
    	double t374 = -t96 * t215 - (-t181 - t208 * t88 / 2 - t211 * t88 / 2 - t90 * t221 / 2 - t87 * t306 * t152 / 2 + t36 * t286 * t69 / 2 - t93 * t206 / 2 - t87 * t316 * t161 / 2 + 0.3e1 / 0.2e1 * t36 * t368 * t120) * t57;
    	double t378 = -2 * t40 * t32;
    	double t397 = -t96 * t250 - (-t145 - t246 * t88 / 2 - t90 * t378 / 2 + t1 * t306 * t152 / 2 - t33 * t286 * t69 / 2 - t93 * t241 / 2 + t1 * t316 * t161 / 2 - 0.3e1 / 0.2e1 * t33 * t368 * t120 + t83) * t57;
    	double t398 = t139 * t139;
    	double t402 = 2 * t32 * t29;
    	double t406 = t20 * t20;
    	double t408 = 0.1e1 / t21 / t406;
    	double t409 = t29 * t408;
    	double t412 = 0.3e1 / 0.4e1 * t300 * t409 * t45;
    	double t413 = t135 * t45;
    	double t432 = 2 * t87 * t29;
    	double t436 = t408 * t45;
    	double t442 = -t139 * t177 - (-t170 * t146 / 2 - t173 * t128 / 2 - t130 * t186 / 2 + 0.3e1 / 0.2e1 * t35 * t402 * t106 - t1 * t402 * t152 / 2 - t136 * t168 / 2 + t35 * t432 * t152 / 2 - 0.3e1 / 0.2e1 * t1 * t432 * t436) * t57;
    	double t451 = 2 * t32 * t47;
    	double t460 = 2 * t87 * t47;
    	double t466 = -t139 * t215 - (-t208 * t128 / 2 - t340 - t130 * t221 / 2 - t87 * t402 * t152 / 2 + t36 * t451 * t69 / 2 - t136 * t206 / 2 - 0.3e1 / 0.2e1 * t87 * t432 * t436 + t36 * t460 * t161 / 2 + t413) * t57;
    	double t490 = -t139 * t250 - (-t181 - t243 * t128 / 2 - t246 * t128 / 2 - t130 * t378 / 2 + t1 * t402 * t152 / 2 - t33 * t451 * t69 / 2 - t136 * t241 / 2 + 0.3e1 / 0.2e1 * t1 * t432 * t436 - t33 * t460 * t161 / 2) * t57;
    	double t491 = t176 * t176;
    	double t495 = 2 * t35 * t29;
    	double t501 = 0.3e1 / 0.4e1 * t252 * t409 * t45;
    	double t517 = 2 * t35 * t47;
    	double t523 = -2 * t1 * t29;
    	double t527 = -2 * t1 * t47;
    	double t533 = -t176 * t215 - (t181 - t208 * t168 / 2 - t48 * t12 * t79 / 2 - t30 * t40 * t282 / 2 - t87 * t495 * t152 / 2 + t36 * t517 * t69 / 2 - t173 * t206 / 2 - 0.3e1 / 0.2e1 * t87 * t523 * t436 + t36 * t527 * t161 / 2) * t57;
    	double t555 = -t176 * t250 - (-t243 * t168 / 2 - t246 * t168 / 2 - t223 + t1 * t495 * t152 / 2 - t33 * t517 * t69 / 2 - t173 * t241 / 2 + 0.3e1 / 0.2e1 * t1 * t523 * t436 - t33 * t527 * t161 / 2 + t413) * t57;
    	double t556 = t214 * t214;
    	double t560 = -2 * t87 * t47;
    	double t593 = -t214 * t250 - (-t243 * t206 / 2 - t246 * t206 / 2 - t208 * t241 / 2 - 3 * t1 * t87 * t29 * t436 - t33 * t560 * t161 / 2 - t211 * t241 / 2 + t1 * t36 * t47 * t161 - 3 * t33 * t36 * t77 * t120) * t57;
    	double t594 = t249 * t249;
    	K[0][0] = -t53 * kb - (-t42 * t58 - t49 * t13 + t68 - t36 * t70 * t69 - t74 + t82 - t83) * t57;
    	K[0][1] = t127;
    	K[0][2] = t167;
    	K[0][3] = t205;
    	K[0][4] = t240;
    	K[0][5] = t277;
    	K[1][0] = t127;
    	K[1][1] = -t278 * kb - (-t90 * t103 - t93 * t88 + t285 + t33 * t286 * t69 - t74 + t293 - t83) * t57;
    	K[1][2] = t325;
    	K[1][3] = t349;
    	K[1][4] = t374;
    	K[1][5] = t397;
    	K[2][0] = t167;
    	K[2][1] = t325;
    	K[2][2] = -t398 * kb - (-t130 * t146 - t136 * t128 + t68 + t87 * t402 * t152 - t74 + t412 - t413) * t57;
    	K[2][3] = t442;
    	K[2][4] = t466;
    	K[2][5] = t490;
    	K[3][0] = t205;
    	K[3][1] = t349;
    	K[3][2] = t442;
    	K[3][3] = -t491 * kb - (-t170 * t186 - t173 * t168 + t285 - t1 * t495 * t152 - t74 + t501 - t413) * t57;
    	K[3][4] = t533;
    	K[3][5] = t555;
    	K[4][0] = t240;
    	K[4][1] = t374;
    	K[4][2] = t466;
    	K[4][3] = t533;
    	K[4][4] = -t556 * kb - (-t208 * t206 - t211 * t206 + t412 + t36 * t560 * t161 - t413 + t82 - t83) * t57;
    	K[4][5] = t593;
    	K[5][0] = t277;
    	K[5][1] = t397;
    	K[5][2] = t490;
    	K[5][3] = t555;
    	K[5][4] = t593;
    	K[5][5] = -t594 * kb - (-t243 * t241 - t246 * t241 + t501 - 2 * t33 * t1 * t47 * t161 - t413 + t293 - t83) * t57;
    }
    
}
