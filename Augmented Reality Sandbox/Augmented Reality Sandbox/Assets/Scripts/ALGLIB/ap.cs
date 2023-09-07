/*************************************************************************
AP library
Copyright (c) 2003-2009 Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is a trial version of the ALGLIB package  licensed
to Bilkent University (Licensee, You).

Only Licensee can use  it  according   to   the  ALGLIB   Trial
License Agreement between Licensor (Sole  Proprietor  Bochkanov
Sergey Anatolyevich) and Licensee.

=============== ALGLIB TRIAL LICENSE AGREEMENT ================

1. Only owners of e-mails below can use this trial version:
* aytek.aman@cs.bilkent.edu.tr
* akaydin@cs.bilkent.edu.tr

2. This trial version can be used only for 30 days  (expiration
date is 19 Feb, 2014).

It can be used only for evaluation  purposes. You can  not  use
it to perform some "real" work.

After this trial period is over, you MUST delete all  copies of
software (including ones made for backup  purposes)  -  or  buy
commercial license for ALGLIB.

3. Following restrictions are applied:
A. Trial version of ALGLIB can not be used by other individuals
   or legal entities. It can not  be  distributed  in  any  way
   (including rent, lease, other means of sharing software).
B. Any application/product  developed  with  trial  version  of
   ALGLIB may only be used for  evaluation  purposes  and  only
   during evaluation period.
   You may not  distribute  applications/products  linked  with
   trial version of ALGLIB.
   After evaluation period  is  over, you should  re-link  such
   applications/products  with  Free  Edition  or  one  of  the
   Commercial editions of ALGLIB - or delete them.
C. You can not remove any  copyright  notice  from  the  Source
   Codes/Binary Files
   
WARRANTIES:

This trial version  is  provided  AS IS,  with  no  warranties,
express or implied.  Licensor  does  NOT  provide  support  and
maintenance for this trial version. If you need warranties  and
support, you should buy commercial version of ALGLIB.

COPYRIGHT:

Title to the ALGLIB and all copies thereof remain with Licensor
The ALGLIB is copyrighted and is protected by Russian copyright
laws and international treaty provisions.  You  will not remove
any copyright notice  from  the  ALGLIB  files.  You  agree  to
prevent any unauthorized  copying  of  the  ALGLIB.  Except  as
expressly provided herein, Licensor does not grant any  express
or implied right to you  under  Licensor  patents,  copyrights,
trademarks, or trade secret information.
>>> END OF LICENSE >>>
*************************************************************************/
using System;
public partial class alglib
{
    /********************************************************************
    Callback definitions for optimizers/fitters/solvers.
    
    Callbacks for unparameterized (general) functions:
    * ndimensional_func         calculates f(arg), stores result to func
    * ndimensional_grad         calculates func = f(arg), 
                                grad[i] = df(arg)/d(arg[i])
    * ndimensional_hess         calculates func = f(arg),
                                grad[i] = df(arg)/d(arg[i]),
                                hess[i,j] = d2f(arg)/(d(arg[i])*d(arg[j]))
    
    Callbacks for systems of functions:
    * ndimensional_fvec         calculates vector function f(arg),
                                stores result to fi
    * ndimensional_jac          calculates f[i] = fi(arg)
                                jac[i,j] = df[i](arg)/d(arg[j])
                                
    Callbacks for  parameterized  functions,  i.e.  for  functions  which 
    depend on two vectors: P and Q.  Gradient  and Hessian are calculated 
    with respect to P only.
    * ndimensional_pfunc        calculates f(p,q),
                                stores result to func
    * ndimensional_pgrad        calculates func = f(p,q),
                                grad[i] = df(p,q)/d(p[i])
    * ndimensional_phess        calculates func = f(p,q),
                                grad[i] = df(p,q)/d(p[i]),
                                hess[i,j] = d2f(p,q)/(d(p[i])*d(p[j]))

    Callbacks for progress reports:
    * ndimensional_rep          reports current position of optimization algo    
    
    Callbacks for ODE solvers:
    * ndimensional_ode_rp       calculates dy/dx for given y[] and x
    
    Callbacks for integrators:
    * integrator1_func          calculates f(x) for given x
                                (additional parameters xminusa and bminusx
                                contain x-a and b-x)
    ********************************************************************/
    public delegate void ndimensional_func (double[] arg, ref double func, object obj);
    public delegate void ndimensional_grad (double[] arg, ref double func, double[] grad, object obj);
    public delegate void ndimensional_hess (double[] arg, ref double func, double[] grad, double[,] hess, object obj);
    
    public delegate void ndimensional_fvec (double[] arg, double[] fi, object obj);
    public delegate void ndimensional_jac  (double[] arg, double[] fi, double[,] jac, object obj);
    
    public delegate void ndimensional_pfunc(double[] p, double[] q, ref double func, object obj);
    public delegate void ndimensional_pgrad(double[] p, double[] q, ref double func, double[] grad, object obj);
    public delegate void ndimensional_phess(double[] p, double[] q, ref double func, double[] grad, double[,] hess, object obj);
    
    public delegate void ndimensional_rep(double[] arg, double func, object obj);

    public delegate void ndimensional_ode_rp (double[] y, double x, double[] dy, object obj);

    public delegate void integrator1_func (double x, double xminusa, double bminusx, ref double f, object obj);

    /********************************************************************
    Class defining a complex number with double precision.
    ********************************************************************/
    public struct complex
    {
        public double x;
        public double y;

        public complex(double _x)
        {
            x = _x;
            y = 0;
        }
        public complex(double _x, double _y)
        {
            x = _x;
            y = _y;
        }
        public static implicit operator complex(double _x)
        {
            return new complex(_x);
        }
        public static bool operator==(complex lhs, complex rhs)
        {
            return ((double)lhs.x==(double)rhs.x) & ((double)lhs.y==(double)rhs.y);
        }
        public static bool operator!=(complex lhs, complex rhs)
        {
            return ((double)lhs.x!=(double)rhs.x) | ((double)lhs.y!=(double)rhs.y);
        }
        public static complex operator+(complex lhs)
        {
            return lhs;
        }
        public static complex operator-(complex lhs)
        {
            return new complex(-lhs.x,-lhs.y);
        }
        public static complex operator+(complex lhs, complex rhs)
        {
            return new complex(lhs.x+rhs.x,lhs.y+rhs.y);
        }
        public static complex operator-(complex lhs, complex rhs)
        {
            return new complex(lhs.x-rhs.x,lhs.y-rhs.y);
        }
        public static complex operator*(complex lhs, complex rhs)
        { 
            return new complex(lhs.x*rhs.x-lhs.y*rhs.y, lhs.x*rhs.y+lhs.y*rhs.x);
        }
        public static complex operator/(complex lhs, complex rhs)
        {
            complex result;
            double e;
            double f;
            if( System.Math.Abs(rhs.y)<System.Math.Abs(rhs.x) )
            {
                e = rhs.y/rhs.x;
                f = rhs.x+rhs.y*e;
                result.x = (lhs.x+lhs.y*e)/f;
                result.y = (lhs.y-lhs.x*e)/f;
            }
            else
            {
                e = rhs.x/rhs.y;
                f = rhs.y+rhs.x*e;
                result.x = (lhs.y+lhs.x*e)/f;
                result.y = (-lhs.x+lhs.y*e)/f;
            }
            return result;
        }
        public override int GetHashCode() 
        { 
            return x.GetHashCode() ^ y.GetHashCode(); 
        }
        public override bool Equals(object obj) 
        { 
            if( obj is byte)
                return Equals(new complex((byte)obj));
            if( obj is sbyte)
                return Equals(new complex((sbyte)obj));
            if( obj is short)
                return Equals(new complex((short)obj));
            if( obj is ushort)
                return Equals(new complex((ushort)obj));
            if( obj is int)
                return Equals(new complex((int)obj));
            if( obj is uint)
                return Equals(new complex((uint)obj));
            if( obj is long)
                return Equals(new complex((long)obj));
            if( obj is ulong)
                return Equals(new complex((ulong)obj));
            if( obj is float)
                return Equals(new complex((float)obj));
            if( obj is double)
                return Equals(new complex((double)obj));
            if( obj is decimal)
                return Equals(new complex((double)(decimal)obj));
            return base.Equals(obj); 
        }    
    }    
    
    /********************************************************************
    Class defining an ALGLIB exception
    ********************************************************************/
    public class alglibexception : System.Exception
    {
        public string msg;
        public alglibexception(string s)
        {
            msg = s;
        }
        
    }
    
    /********************************************************************
    ALGLIB object, parent  class  for  all  internal  AlgoPascal  objects
    managed by ALGLIB.
    
    Any internal AlgoPascal object inherits from this class.
    
    User-visible objects inherit from alglibobject (see below).
    ********************************************************************/
    public abstract class apobject
    {
        public abstract void init();
        public abstract apobject make_copy();
    }
    
    /********************************************************************
    ALGLIB object, parent class for all user-visible objects  managed  by
    ALGLIB.
    
    Methods:
        _deallocate()       deallocation:
                            * in managed ALGLIB it does nothing
                            * in native ALGLIB it clears  dynamic  memory
                              being  hold  by  object  and  sets internal
                              reference to null.
    ********************************************************************/
    public abstract class alglibobject
    {
        public virtual void _deallocate() {}
    }
    
    /********************************************************************
    Deallocation of ALGLIB object:
    * in managed ALGLIB this method just sets refence to null
    * in native ALGLIB call of this method:
      1) clears dynamic memory being hold by  object  and  sets  internal
         reference to null.
      2) sets to null variable being passed to this method
    
    IMPORTANT (1): in  native  edition  of  ALGLIB,  obj becomes unusable
                   after this call!!!  It  is  possible  to  save  a copy
                   of reference in another variable (original variable is
                   set to null), but any attempt to work with this object
                   will crash your program.
    
    IMPORTANT (2): memory ownen by object will be recycled by GC  in  any
                   case. This method just enforced IMMEDIATE deallocation.
    ********************************************************************/
    public static void deallocateimmediately<T>(ref T obj) where T : alglib.alglibobject
    {
        obj._deallocate();
        obj = null;
    }

    /********************************************************************
    Allocation counter:
    * in managed ALGLIB it always returns 0 (dummy code)
    * in native ALGLIB it returns current value of the allocation counter
      (if it was activated)
    ********************************************************************/
    public static long alloc_counter()
    {
        return 0;
    }
    
    /********************************************************************
    Activization of the allocation counter:
    * in managed ALGLIB it does nothing (dummy code)
    * in native ALGLIB it turns on allocation counting.
    ********************************************************************/
    public static void alloc_counter_activate()
    {
    }
    
    /********************************************************************
    reverse communication structure
    ********************************************************************/
    public class rcommstate : apobject
    {
        public rcommstate()
        {
            init();
        }
        public override void init()
        {
            stage = -1;
            ia = new int[0];
            ba = new bool[0];
            ra = new double[0];
            ca = new alglib.complex[0];
        }
        public override apobject make_copy()
        {
            rcommstate result = new rcommstate();
            result.stage = stage;
            result.ia = (int[])ia.Clone();
            result.ba = (bool[])ba.Clone();
            result.ra = (double[])ra.Clone();
            result.ca = (alglib.complex[])ca.Clone();
            return result;
        }
        public int stage;
        public int[] ia;
        public bool[] ba;
        public double[] ra;
        public alglib.complex[] ca;
    };

    /********************************************************************
    internal functions
    ********************************************************************/
    public class ap
    {
        public static int len<T>(T[] a)
        { return a.Length; }
        public static int rows<T>(T[,] a)
        { return a.GetLength(0); }
        public static int cols<T>(T[,] a)
        { return a.GetLength(1); }
        public static void swap<T>(ref T a, ref T b)
        {
            T t = a;
            a = b;
            b = t;
        }
        
        public static void assert(bool cond, string s)
        {
            if( !cond )
                throw new alglibexception(s);
        }
        
        public static void assert(bool cond)
        {
            assert(cond, "ALGLIB: assertion failed");
        }
        
        /****************************************************************
        returns dps (digits-of-precision) value corresponding to threshold.
        dps(0.9)  = dps(0.5)  = dps(0.1) = 0
        dps(0.09) = dps(0.05) = dps(0.01) = 1
        and so on
        ****************************************************************/
        public static int threshold2dps(double threshold)
        {
            int result = 0;
            double t;
            for (result = 0, t = 1; t / 10 > threshold*(1+1E-10); result++, t /= 10) ;
            return result;
        }

        /****************************************************************
        prints formatted complex
        ****************************************************************/
        public static string format(complex a, int _dps)
        {
            int dps = Math.Abs(_dps);
            string fmt = _dps>=0 ? "F" : "E";
            string fmtx = String.Format("{{0:"+fmt+"{0}}}", dps);
            string fmty = String.Format("{{0:"+fmt+"{0}}}", dps);
            string result = String.Format(fmtx, a.x) + (a.y >= 0 ? "+" : "-") + String.Format(fmty, Math.Abs(a.y)) + "i";
            result = result.Replace(',', '.');
            return result;
        }

        /****************************************************************
        prints formatted array
        ****************************************************************/
        public static string format(bool[] a)
        {
            string[] result = new string[len(a)];
            int i;
            for(i=0; i<len(a); i++)
                if( a[i] )
                    result[i] = "true";
                else
                    result[i] = "false";
            return "{"+String.Join(",",result)+"}";
        }
        
        /****************************************************************
        prints formatted array
        ****************************************************************/
        public static string format(int[] a)
        {
            string[] result = new string[len(a)];
            int i;
            for (i = 0; i < len(a); i++)
                result[i] = a[i].ToString();
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted array
        ****************************************************************/
        public static string format(double[] a, int _dps)
        {
            int dps = Math.Abs(_dps);
            string sfmt = _dps >= 0 ? "F" : "E";
            string fmt = String.Format("{{0:" + sfmt + "{0}}}", dps);
            string[] result = new string[len(a)];
            int i;
            for (i = 0; i < len(a); i++)
            {
                result[i] = String.Format(fmt, a[i]);
                result[i] = result[i].Replace(',', '.');
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted array
        ****************************************************************/
        public static string format(complex[] a, int _dps)
        {
            int dps = Math.Abs(_dps);
            string fmt = _dps >= 0 ? "F" : "E";
            string fmtx = String.Format("{{0:"+fmt+"{0}}}", dps);
            string fmty = String.Format("{{0:"+fmt+"{0}}}", dps);
            string[] result = new string[len(a)];
            int i;
            for (i = 0; i < len(a); i++)
            {
                result[i] = String.Format(fmtx, a[i].x) + (a[i].y >= 0 ? "+" : "-") + String.Format(fmty, Math.Abs(a[i].y)) + "i";
                result[i] = result[i].Replace(',', '.');
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted matrix
        ****************************************************************/
        public static string format(bool[,] a)
        {
            int i, j, m, n;
            n = cols(a);
            m = rows(a);
            bool[] line = new bool[n];
            string[] result = new string[m];
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    line[j] = a[i, j];
                result[i] = format(line);
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted matrix
        ****************************************************************/
        public static string format(int[,] a)
        {
            int i, j, m, n;
            n = cols(a);
            m = rows(a);
            int[] line = new int[n];
            string[] result = new string[m];
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    line[j] = a[i, j];
                result[i] = format(line);
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted matrix
        ****************************************************************/
        public static string format(double[,] a, int dps)
        {
            int i, j, m, n;
            n = cols(a);
            m = rows(a);
            double[] line = new double[n];
            string[] result = new string[m];
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    line[j] = a[i, j];
                result[i] = format(line, dps);
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        prints formatted matrix
        ****************************************************************/
        public static string format(complex[,] a, int dps)
        {
            int i, j, m, n;
            n = cols(a);
            m = rows(a);
            complex[] line = new complex[n];
            string[] result = new string[m];
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    line[j] = a[i, j];
                result[i] = format(line, dps);
            }
            return "{" + String.Join(",", result) + "}";
        }

        /****************************************************************
        checks that matrix is symmetric.
        max|A-A^T| is calculated; if it is within 1.0E-14 of max|A|,
        matrix is considered symmetric
        ****************************************************************/
        public static bool issymmetric(double[,] a)
        {
            int i, j, n;
            double err, mx, v1, v2;
            if( rows(a)!=cols(a) )
                return false;
            n = rows(a);
            if( n==0 )
                return true;
            mx = 0;
            err = 0;
            for( i=0; i<n; i++)
            {
                for(j=i+1; j<n; j++)
                {
                    v1 = a[i,j];
                    v2 = a[j,i];
                    if( !math.isfinite(v1) )
                        return false;
                    if( !math.isfinite(v2) )
                        return false;
                    err = Math.Max(err, Math.Abs(v1-v2));
                    mx  = Math.Max(mx,  Math.Abs(v1));
                    mx  = Math.Max(mx,  Math.Abs(v2));
                }
                v1 = a[i,i];
                if( !math.isfinite(v1) )
                    return false;
                mx = Math.Max(mx, Math.Abs(v1));
            }
            if( mx==0 )
                return true;
            return err/mx<=1.0E-14;
        }
        
        /****************************************************************
        checks that matrix is Hermitian.
        max|A-A^H| is calculated; if it is within 1.0E-14 of max|A|,
        matrix is considered Hermitian
        ****************************************************************/
        public static bool ishermitian(complex[,] a)
        {
            int i, j, n;
            double err, mx;
            complex v1, v2, vt;
            if( rows(a)!=cols(a) )
                return false;
            n = rows(a);
            if( n==0 )
                return true;
            mx = 0;
            err = 0;
            for( i=0; i<n; i++)
            {
                for(j=i+1; j<n; j++)
                {
                    v1 = a[i,j];
                    v2 = a[j,i];
                    if( !math.isfinite(v1.x) )
                        return false;
                    if( !math.isfinite(v1.y) )
                        return false;
                    if( !math.isfinite(v2.x) )
                        return false;
                    if( !math.isfinite(v2.y) )
                        return false;
                    vt.x = v1.x-v2.x;
                    vt.y = v1.y+v2.y;
                    err = Math.Max(err, math.abscomplex(vt));
                    mx  = Math.Max(mx,  math.abscomplex(v1));
                    mx  = Math.Max(mx,  math.abscomplex(v2));
                }
                v1 = a[i,i];
                if( !math.isfinite(v1.x) )
                    return false;
                if( !math.isfinite(v1.y) )
                    return false;
                err = Math.Max(err, Math.Abs(v1.y));
                mx = Math.Max(mx, math.abscomplex(v1));
            }
            if( mx==0 )
                return true;
            return err/mx<=1.0E-14;
        }
        
        
        /****************************************************************
        Forces symmetricity by copying upper half of A to the lower one
        ****************************************************************/
        public static bool forcesymmetric(double[,] a)
        {
            int i, j, n;
            if( rows(a)!=cols(a) )
                return false;
            n = rows(a);
            if( n==0 )
                return true;
            for( i=0; i<n; i++)
                for(j=i+1; j<n; j++)
                    a[i,j] = a[j,i];
            return true;
        }
        
        /****************************************************************
        Forces Hermiticity by copying upper half of A to the lower one
        ****************************************************************/
        public static bool forcehermitian(complex[,] a)
        {
            int i, j, n;
            complex v;
            if( rows(a)!=cols(a) )
                return false;
            n = rows(a);
            if( n==0 )
                return true;
            for( i=0; i<n; i++)
                for(j=i+1; j<n; j++)
                {
                    v = a[j,i];
                    a[i,j].x = v.x;
                    a[i,j].y = -v.y;
                }
            return true;
        }
    };
    
    /********************************************************************
    math functions
    ********************************************************************/
    public class math
    {
        //public static System.Random RndObject = new System.Random(System.DateTime.Now.Millisecond);
        public static System.Random rndobject = new System.Random(System.DateTime.Now.Millisecond + 1000*System.DateTime.Now.Second + 60*1000*System.DateTime.Now.Minute);

        public const double machineepsilon = 5E-16;
        public const double maxrealnumber = 1E300;
        public const double minrealnumber = 1E-300;
        
        public static bool isfinite(double d)
        {
            return !System.Double.IsNaN(d) && !System.Double.IsInfinity(d);
        }
        
        public static double randomreal()
        {
            double r = 0;
            lock(rndobject){ r = rndobject.NextDouble(); }
            return r;
        }
        public static int randominteger(int N)
        {
            int r = 0;
            lock(rndobject){ r = rndobject.Next(N); }
            return r;
        }
        public static double sqr(double X)
        {
            return X*X;
        }        
        public static double abscomplex(complex z)
        {
            double w;
            double xabs;
            double yabs;
            double v;
    
            xabs = System.Math.Abs(z.x);
            yabs = System.Math.Abs(z.y);
            w = xabs>yabs ? xabs : yabs;
            v = xabs<yabs ? xabs : yabs; 
            if( v==0 )
                return w;
            else
            {
                double t = v/w;
                return w*System.Math.Sqrt(1+t*t);
            }
        }
        public static complex conj(complex z)
        {
            return new complex(z.x, -z.y); 
        }    
        public static complex csqr(complex z)
        {
            return new complex(z.x*z.x-z.y*z.y, 2*z.x*z.y); 
        }

    }
    
    
    /********************************************************************
    serializer object (should not be used directly)
    ********************************************************************/
    public class serializer
    {
        enum SMODE { DEFAULT, ALLOC, TO_STRING, FROM_STRING };
        private const int SER_ENTRIES_PER_ROW = 5;
        private const int SER_ENTRY_LENGTH    = 11;
        
        private SMODE mode;
        private int entries_needed;
        private int entries_saved;
        private int bytes_asked;
        private int bytes_written;
        private int bytes_read;
        private char[] out_str;
        private char[] in_str;
        
        public serializer()
        {
            mode = SMODE.DEFAULT;
            entries_needed = 0;
            bytes_asked = 0;
        }

        public void alloc_start()
        {
            entries_needed = 0;
            bytes_asked = 0;
            mode = SMODE.ALLOC;
        }

        public void alloc_entry()
        {
            if( mode!=SMODE.ALLOC )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            entries_needed++;
        }

        private int get_alloc_size()
        {
            int rows, lastrowsize, result;
            
            // check and change mode
            if( mode!=SMODE.ALLOC )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            
            // if no entries needes (degenerate case)
            if( entries_needed==0 )
            {
                bytes_asked = 1;
                return bytes_asked;
            }
            
            // non-degenerate case
            rows = entries_needed/SER_ENTRIES_PER_ROW;
            lastrowsize = SER_ENTRIES_PER_ROW;
            if( entries_needed%SER_ENTRIES_PER_ROW!=0 )
            {
                lastrowsize = entries_needed%SER_ENTRIES_PER_ROW;
                rows++;
            }
            
            // calculate result size
            result  = ((rows-1)*SER_ENTRIES_PER_ROW+lastrowsize)*SER_ENTRY_LENGTH;
            result +=  (rows-1)*(SER_ENTRIES_PER_ROW-1)+(lastrowsize-1);
            result += rows*2;
            bytes_asked = result;
            return result;
        }

        public void sstart_str()
        {
            int allocsize = get_alloc_size();
            
            // check and change mode
            if( mode!=SMODE.ALLOC )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            mode = SMODE.TO_STRING;
            
            // other preparations
            out_str = new char[allocsize];
            entries_saved = 0;
            bytes_written = 0;
        }

        public void ustart_str(string s)
        {
            // check and change mode
            if( mode!=SMODE.DEFAULT )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            mode = SMODE.FROM_STRING;
            
            in_str = s.ToCharArray();
            bytes_read = 0;
        }

        public void serialize_bool(bool v)
        {
            if( mode!=SMODE.TO_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            bool2str(v, out_str, ref bytes_written);
            entries_saved++;
            if( entries_saved%SER_ENTRIES_PER_ROW!=0 )
            {
                out_str[bytes_written] = ' ';
                bytes_written++;
            }
            else
            {
                out_str[bytes_written+0] = '\r';
                out_str[bytes_written+1] = '\n';
                bytes_written+=2;
            }            
        }

        public void serialize_int(int v)
        {
            if( mode!=SMODE.TO_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            int2str(v, out_str, ref bytes_written);
            entries_saved++;
            if( entries_saved%SER_ENTRIES_PER_ROW!=0 )
            {
                out_str[bytes_written] = ' ';
                bytes_written++;
            }
            else
            {
                out_str[bytes_written+0] = '\r';
                out_str[bytes_written+1] = '\n';
                bytes_written+=2;
            }
        }

        public void serialize_double(double v)
        {
            if( mode!=SMODE.TO_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            double2str(v, out_str, ref bytes_written);
            entries_saved++;
            if( entries_saved%SER_ENTRIES_PER_ROW!=0 )
            {
                out_str[bytes_written] = ' ';
                bytes_written++;
            }
            else
            {
                out_str[bytes_written+0] = '\r';
                out_str[bytes_written+1] = '\n';
                bytes_written+=2;
            }
        }

        public bool unserialize_bool()
        {
            if( mode!=SMODE.FROM_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            return str2bool(in_str, ref bytes_read);
        }

        public int unserialize_int()
        {
            if( mode!=SMODE.FROM_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            return str2int(in_str, ref bytes_read);
        }

        public double unserialize_double()
        {
            if( mode!=SMODE.FROM_STRING )
                throw new alglib.alglibexception("ALGLIB: internal error during (un)serialization");
            return str2double(in_str, ref bytes_read);
        }

        public void stop()
        {
        }

        public string get_string()
        {
            return new string(out_str, 0, bytes_written);
        }


        /************************************************************************
        This function converts six-bit value (from 0 to 63)  to  character  (only
        digits, lowercase and uppercase letters, minus and underscore are used).

        If v is negative or greater than 63, this function returns '?'.
        ************************************************************************/
        private static char[] _sixbits2char_tbl = new char[64]{ 
                '0', '1', '2', '3', '4', '5', '6', '7',
                '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
                'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
                'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
                'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 
                'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 
                'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 
                'u', 'v', 'w', 'x', 'y', 'z', '-', '_' };
        private static char sixbits2char(int v)
        {
            if( v<0 || v>63 )
                return '?';
            return _sixbits2char_tbl[v];
        }
        
        /************************************************************************
        This function converts character to six-bit value (from 0 to 63).

        This function is inverse of ae_sixbits2char()
        If c is not correct character, this function returns -1.
        ************************************************************************/
        private static int[] _char2sixbits_tbl = new int[128] {
            -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, 62, -1, -1,
             0,  1,  2,  3,  4,  5,  6,  7,
             8,  9, -1, -1, -1, -1, -1, -1,
            -1, 10, 11, 12, 13, 14, 15, 16,
            17, 18, 19, 20, 21, 22, 23, 24,
            25, 26, 27, 28, 29, 30, 31, 32,
            33, 34, 35, -1, -1, -1, -1, 63,
            -1, 36, 37, 38, 39, 40, 41, 42,
            43, 44, 45, 46, 47, 48, 49, 50,
            51, 52, 53, 54, 55, 56, 57, 58,
            59, 60, 61, -1, -1, -1, -1, -1 };
        private static int char2sixbits(char c)
        {
            return (c>=0 && c<127) ? _char2sixbits_tbl[c] : -1;
        }
        
        /************************************************************************
        This function converts three bytes (24 bits) to four six-bit values 
        (24 bits again).

        src         array
        src_offs    offset of three-bytes chunk
        dst         array for ints
        dst_offs    offset of four-ints chunk
        ************************************************************************/
        private static void threebytes2foursixbits(byte[] src, int src_offs, int[] dst, int dst_offs)
        {
            dst[dst_offs+0] =  src[src_offs+0] & 0x3F;
            dst[dst_offs+1] = (src[src_offs+0]>>6) | ((src[src_offs+1]&0x0F)<<2);
            dst[dst_offs+2] = (src[src_offs+1]>>4) | ((src[src_offs+2]&0x03)<<4);
            dst[dst_offs+3] =  src[src_offs+2]>>2;
        }

        /************************************************************************
        This function converts four six-bit values (24 bits) to three bytes
        (24 bits again).

        src         pointer to four ints
        src_offs    offset of the chunk
        dst         pointer to three bytes
        dst_offs    offset of the chunk
        ************************************************************************/
        private static void foursixbits2threebytes(int[] src, int src_offs, byte[] dst, int dst_offs)
        {
            dst[dst_offs+0] =      (byte)(src[src_offs+0] | ((src[src_offs+1]&0x03)<<6));
            dst[dst_offs+1] = (byte)((src[src_offs+1]>>2) | ((src[src_offs+2]&0x0F)<<4));
            dst[dst_offs+2] = (byte)((src[src_offs+2]>>4) |  (src[src_offs+3]<<2));
        }

        /************************************************************************
        This function serializes boolean value into buffer

        v           boolean value to be serialized
        buf         buffer, at least 11 characters wide
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.
        ************************************************************************/
        private static void bool2str(bool v, char[] buf, ref int offs)
        {
            char c = v ? '1' : '0';
            int i;
            for(i=0; i<SER_ENTRY_LENGTH; i++)
                buf[offs+i] = c;
            offs += SER_ENTRY_LENGTH;
        }

        /************************************************************************
        This function unserializes boolean value from buffer

        buf         buffer which contains value; leading spaces/tabs/newlines are 
                    ignored, traling spaces/tabs/newlines are treated as  end  of
                    the boolean value.
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.

        This function raises an error in case unexpected symbol is found
        ************************************************************************/
        private static bool str2bool(char[] buf, ref int offs)
        {
            bool was0, was1;
            string emsg = "ALGLIB: unable to read boolean value from stream";
            
            was0 = false;
            was1 = false;
            while( buf[offs]==' ' || buf[offs]=='\t' || buf[offs]=='\n' || buf[offs]=='\r' )
                offs++;
            while( buf[offs]!=' ' && buf[offs]!='\t' && buf[offs]!='\n' && buf[offs]!='\r' && buf[offs]!=0 )
            {
                if( buf[offs]=='0' )
                {
                    was0 = true;
                    offs++;
                    continue;
                }
                if( buf[offs]=='1' )
                {
                    was1 = true;
                    offs++;
                    continue;
                }
                throw new alglib.alglibexception(emsg);
            }
            if( (!was0) && (!was1) )
                throw new alglib.alglibexception(emsg);
            if( was0 && was1 )
                throw new alglib.alglibexception(emsg);
            return was1 ? true : false;
        }

        /************************************************************************
        This function serializes integer value into buffer

        v           integer value to be serialized
        buf         buffer, at least 11 characters wide 
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.

        This function raises an error in case unexpected symbol is found
        ************************************************************************/
        private static void int2str(int v, char[] buf, ref int offs)
        {
            int i;
            byte[] _bytes = System.BitConverter.GetBytes((int)v);
            byte[]  bytes = new byte[9];
            int[] sixbits = new int[12];
            byte c;
            
            //
            // copy v to array of bytes, sign extending it and 
            // converting to little endian order. Additionally, 
            // we set 9th byte to zero in order to simplify 
            // conversion to six-bit representation
            //
            if( !System.BitConverter.IsLittleEndian )
                System.Array.Reverse(_bytes);
            c = v<0 ? (byte)0xFF : (byte)0x00;
            for(i=0; i<sizeof(int); i++)
                bytes[i] = _bytes[i];
            for(i=sizeof(int); i<8; i++)
                bytes[i] = c;
            bytes[8] = 0;
            
            //
            // convert to six-bit representation, output
            //
            // NOTE: last 12th element of sixbits is always zero, we do not output it
            //
            threebytes2foursixbits(bytes, 0, sixbits, 0);
            threebytes2foursixbits(bytes, 3, sixbits, 4);
            threebytes2foursixbits(bytes, 6, sixbits, 8);        
            for(i=0; i<SER_ENTRY_LENGTH; i++)
                buf[offs+i] = sixbits2char(sixbits[i]);
            offs += SER_ENTRY_LENGTH;
        }

        /************************************************************************
        This function unserializes integer value from string

        buf         buffer which contains value; leading spaces/tabs/newlines are 
                    ignored, traling spaces/tabs/newlines are treated as  end  of
                    the integer value.
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.

        This function raises an error in case unexpected symbol is found
        ************************************************************************/
        private static int str2int(char[] buf, ref int offs)
        {
            string emsg =       "ALGLIB: unable to read integer value from stream";
            string emsg3264 =   "ALGLIB: unable to read integer value from stream (value does not fit into 32 bits)";
            int[] sixbits = new int[12];
            byte[] bytes = new byte[9];
            byte[] _bytes = new byte[sizeof(int)];
            int sixbitsread, i;
            byte c;
            
            // 
            // 1. skip leading spaces
            // 2. read and decode six-bit digits
            // 3. set trailing digits to zeros
            // 4. convert to little endian 64-bit integer representation
            // 5. check that we fit into int
            // 6. convert to big endian representation, if needed
            //
            sixbitsread = 0;
            while( buf[offs]==' ' || buf[offs]=='\t' || buf[offs]=='\n' || buf[offs]=='\r' )
                offs++;
            while( buf[offs]!=' ' && buf[offs]!='\t' && buf[offs]!='\n' && buf[offs]!='\r' && buf[offs]!=0 )
            {
                int d;
                d = char2sixbits(buf[offs]);
                if( d<0 || sixbitsread>=SER_ENTRY_LENGTH )
                    throw new alglib.alglibexception(emsg);
                sixbits[sixbitsread] = d;
                sixbitsread++;
                offs++;
            }
            if( sixbitsread==0 )
                throw new alglib.alglibexception(emsg);
            for(i=sixbitsread; i<12; i++)
                sixbits[i] = 0;
            foursixbits2threebytes(sixbits, 0, bytes, 0);
            foursixbits2threebytes(sixbits, 4, bytes, 3);
            foursixbits2threebytes(sixbits, 8, bytes, 6);
            c = (bytes[sizeof(int)-1] & 0x80)!=0 ? (byte)0xFF : (byte)0x00;
            for(i=sizeof(int); i<8; i++)
                if( bytes[i]!=c )
                    throw new alglib.alglibexception(emsg3264);
            for(i=0; i<sizeof(int); i++)
                _bytes[i] = bytes[i];        
            if( !System.BitConverter.IsLittleEndian )
                System.Array.Reverse(_bytes);
            return System.BitConverter.ToInt32(_bytes,0);
        }    
        
        
        /************************************************************************
        This function serializes double value into buffer

        v           double value to be serialized
        buf         buffer, at least 11 characters wide 
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.
        ************************************************************************/
        private static void double2str(double v, char[] buf, ref int offs)
        {
            int i;
            int[] sixbits = new int[12];
            byte[] bytes = new byte[9];

            //
            // handle special quantities
            //
            if( System.Double.IsNaN(v) )
            {
                buf[offs+0] = '.';
                buf[offs+1] = 'n';
                buf[offs+2] = 'a';
                buf[offs+3] = 'n';
                buf[offs+4] = '_';
                buf[offs+5] = '_';
                buf[offs+6] = '_';
                buf[offs+7] = '_';
                buf[offs+8] = '_';
                buf[offs+9] = '_';
                buf[offs+10] = '_';
                offs += SER_ENTRY_LENGTH;
                return;
            }
            if( System.Double.IsPositiveInfinity(v) )
            {
                buf[offs+0] = '.';
                buf[offs+1] = 'p';
                buf[offs+2] = 'o';
                buf[offs+3] = 's';
                buf[offs+4] = 'i';
                buf[offs+5] = 'n';
                buf[offs+6] = 'f';
                buf[offs+7] = '_';
                buf[offs+8] = '_';
                buf[offs+9] = '_';
                buf[offs+10] = '_';
                offs += SER_ENTRY_LENGTH;
                return;
            }
            if( System.Double.IsNegativeInfinity(v) )
            {
                buf[offs+0] = '.';
                buf[offs+1] = 'n';
                buf[offs+2] = 'e';
                buf[offs+3] = 'g';
                buf[offs+4] = 'i';
                buf[offs+5] = 'n';
                buf[offs+6] = 'f';
                buf[offs+7] = '_';
                buf[offs+8] = '_';
                buf[offs+9] = '_';
                buf[offs+10] = '_';
                offs += SER_ENTRY_LENGTH;
                return;
            }
            
            //
            // process general case:
            // 1. copy v to array of chars
            // 2. set 9th byte to zero in order to simplify conversion to six-bit representation
            // 3. convert to little endian (if needed)
            // 4. convert to six-bit representation
            //    (last 12th element of sixbits is always zero, we do not output it)
            //
            byte[] _bytes = System.BitConverter.GetBytes((double)v);
            if( !System.BitConverter.IsLittleEndian )
                System.Array.Reverse(_bytes);
            for(i=0; i<sizeof(double); i++)
                bytes[i] = _bytes[i];
            for(i=sizeof(double); i<9; i++)
                bytes[i] = 0;
            threebytes2foursixbits(bytes, 0, sixbits, 0);
            threebytes2foursixbits(bytes, 3, sixbits, 4);
            threebytes2foursixbits(bytes, 6, sixbits, 8);
            for(i=0; i<SER_ENTRY_LENGTH; i++)
                buf[offs+i] = sixbits2char(sixbits[i]);
            offs += SER_ENTRY_LENGTH;
        }

        /************************************************************************
        This function unserializes double value from string

        buf         buffer which contains value; leading spaces/tabs/newlines are 
                    ignored, traling spaces/tabs/newlines are treated as  end  of
                    the double value.
        offs        offset in the buffer
        
        after return from this function, offs points to the char's past the value
        being read.

        This function raises an error in case unexpected symbol is found
        ************************************************************************/
        private static double str2double(char[] buf, ref int offs)
        {
            string emsg = "ALGLIB: unable to read double value from stream";
            int[] sixbits = new int[12];
            byte[]  bytes = new byte[9];
            byte[] _bytes = new byte[sizeof(double)];
            int sixbitsread, i;
            
            
            // 
            // skip leading spaces
            //
            while( buf[offs]==' ' || buf[offs]=='\t' || buf[offs]=='\n' || buf[offs]=='\r' )
                offs++;
            
              
            //
            // Handle special cases
            //
            if( buf[offs]=='.' )
            {
                string s = new string(buf, offs, SER_ENTRY_LENGTH);
                if( s==".nan_______" )
                {
                    offs += SER_ENTRY_LENGTH;
                    return System.Double.NaN;
                }
                if( s==".posinf____" )
                {
                    offs += SER_ENTRY_LENGTH;
                    return System.Double.PositiveInfinity;
                }
                if( s==".neginf____" )
                {
                    offs += SER_ENTRY_LENGTH;
                    return System.Double.NegativeInfinity;
                }
                throw new alglib.alglibexception(emsg);
            }
            
            // 
            // General case:
            // 1. read and decode six-bit digits
            // 2. check that all 11 digits were read
            // 3. set last 12th digit to zero (needed for simplicity of conversion)
            // 4. convert to 8 bytes
            // 5. convert to big endian representation, if needed
            //
            sixbitsread = 0;
            while( buf[offs]!=' ' && buf[offs]!='\t' && buf[offs]!='\n' && buf[offs]!='\r' && buf[offs]!=0 )
            {
                int d;
                d = char2sixbits(buf[offs]);
                if( d<0 || sixbitsread>=SER_ENTRY_LENGTH )
                    throw new alglib.alglibexception(emsg);
                sixbits[sixbitsread] = d;
                sixbitsread++;
                offs++;
            }
            if( sixbitsread!=SER_ENTRY_LENGTH )
                throw new alglib.alglibexception(emsg);
            sixbits[SER_ENTRY_LENGTH] = 0;
            foursixbits2threebytes(sixbits, 0, bytes, 0);
            foursixbits2threebytes(sixbits, 4, bytes, 3);
            foursixbits2threebytes(sixbits, 8, bytes, 6);
            for(i=0; i<sizeof(double); i++)
                _bytes[i] = bytes[i];        
            if( !System.BitConverter.IsLittleEndian )
                System.Array.Reverse(_bytes);        
            return System.BitConverter.ToDouble(_bytes,0);
        }
    }

    /*
     * Parts of alglib.smp class which are shared with GPL version of ALGLIB
     */
    public partial class smp
    {
        #pragma warning disable 420
        public const int AE_LOCK_CYCLES = 512;
        public const int AE_LOCK_TESTS_BEFORE_YIELD = 16;
        
        /*
         * This variable is used to perform spin-wait loops in a platform-independent manner
         * (loops which should work same way on Mono and Microsoft NET). You SHOULD NEVER
         * change this field - it must be zero during all program life.
         */
        public static volatile int never_change_it = 0;
        
        /*************************************************************************
        Lock.

        This class provides lightweight spin lock
        *************************************************************************/
        public class ae_lock
        {
            public volatile int is_locked;
        }

        /********************************************************************
        Shared pool: data structure used to provide thread-safe access to pool
        of temporary variables.
        ********************************************************************/
        public class sharedpoolentry
        {
            public apobject obj;
            public sharedpoolentry next_entry;
        }
        public class shared_pool : apobject
        {
            /* lock object which protects pool */
            public ae_lock pool_lock;
    
            /* seed object (used to create new instances of temporaries) */
            public volatile apobject seed_object;
            
            /*
             * list of recycled OBJECTS:
             * 1. entries in this list store pointers to recycled objects
             * 2. every time we retrieve object, we retrieve first entry from this list,
             *    move it to recycled_entries and return its obj field to caller/
             */
            public volatile sharedpoolentry recycled_objects;
            
            /* 
             * list of recycled ENTRIES:
             * 1. this list holds entries which are not used to store recycled objects;
             *    every time recycled object is retrieved, its entry is moved to this list.
             * 2. every time object is recycled, we try to fetch entry for him from this list
             *    before allocating it with malloc()
             */
            public volatile sharedpoolentry recycled_entries;
            
            /* enumeration pointer, points to current recycled object*/
            public volatile sharedpoolentry enumeration_counter;
            
            /* constructor */
            public shared_pool()
            {
                ae_init_lock(ref pool_lock);
            }
            
            /* initializer - creation of empty pool */
            public override void init()
            {
                seed_object = null;
                recycled_objects = null;
                recycled_entries = null;
                enumeration_counter = null;
            }
            
            /* copy constructor (it is NOT thread-safe) */
            public override apobject make_copy()
            {
                sharedpoolentry ptr, buf;
                shared_pool result = new shared_pool();
                
                /* create lock */
                ae_init_lock(ref result.pool_lock);
    
                /* copy seed object */
                if( seed_object!=null )
                    result.seed_object = seed_object.make_copy();
                
                /*
                 * copy recycled objects:
                 * 1. copy to temporary list (objects are inserted to beginning, order is reversed)
                 * 2. copy temporary list to output list (order is restored back to normal)
                 */
                buf = null;
                for(ptr=recycled_objects; ptr!=null; ptr=ptr.next_entry)
                {
                    sharedpoolentry tmp = new sharedpoolentry();
                    tmp.obj =  ptr.obj.make_copy();
                    tmp.next_entry = buf;
                    buf = tmp;
                }
                result.recycled_objects = null;
                for(ptr=buf; ptr!=null;)
                {
                    sharedpoolentry next_ptr = ptr.next_entry;
                    ptr.next_entry = result.recycled_objects;
                    result.recycled_objects = ptr;
                    ptr = next_ptr;
                }
    
                /* recycled entries are not copied because they do not store any information */
                result.recycled_entries = null;
    
                /* enumeration counter is reset on copying */
                result.enumeration_counter = null;
    
                return result;
            }
        }
        

        /************************************************************************
        This function performs given number of spin-wait iterations
        ************************************************************************/
        public static void ae_spin_wait(int cnt)
        {
            /*
             * these strange operations with ae_never_change_it are necessary to
             * prevent compiler optimization of the loop.
             */
            int i;
            
            /* very unlikely because no one will wait for such amount of cycles */
            if( cnt>0x12345678 )
                never_change_it = cnt%10;
            
            /* spin wait, test condition which will never be true */
            for(i=0; i<cnt; i++)
                if( never_change_it>0 )
                    never_change_it--;
        }


        /************************************************************************
        This function causes the calling thread to relinquish the CPU. The thread
        is moved to the end of the queue and some other thread gets to run.
        ************************************************************************/
        public static void ae_yield()
        {
            System.Threading.Thread.Sleep(0);
        }

        /************************************************************************
        This function initializes ae_lock structure and sets lock in a free mode.
        ************************************************************************/
        public static void ae_init_lock(ref ae_lock obj)
        {
            obj = new ae_lock();
            obj.is_locked = 0;
        }


        /************************************************************************
        This function acquires lock. In case lock is busy, we perform several
        iterations inside tight loop before trying again.
        ************************************************************************/
        public static void ae_acquire_lock(ae_lock obj)
        {
            int cnt = 0;
            for(;;)
            {
                if( System.Threading.Interlocked.CompareExchange(ref obj.is_locked, 1, 0)==0 )
                    return;
                ae_spin_wait(AE_LOCK_CYCLES);
                cnt++;
                if( cnt%AE_LOCK_TESTS_BEFORE_YIELD==0 )
                    ae_yield();
            }
        }


        /************************************************************************
        This function releases lock.
        ************************************************************************/
        public static void ae_release_lock(ae_lock obj)
        {
            System.Threading.Interlocked.Exchange(ref obj.is_locked, 0);
        }


        /************************************************************************
        This function frees ae_lock structure.
        ************************************************************************/
        public static void ae_free_lock(ref ae_lock obj)
        {
            obj = null;
        }
        
        
        /************************************************************************
        This function returns True, if internal seed object was set.  It  returns
        False for un-seeded pool.

        dst                 destination pool (initialized by constructor function)

        NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
              you should NOT call it when lock can be used by another thread.
        ************************************************************************/
        public static bool ae_shared_pool_is_initialized(shared_pool dst)
        {
            return dst.seed_object!=null;
        }


        /************************************************************************
        This function sets internal seed object. All objects owned by the pool
        (current seed object, recycled objects) are automatically freed.

        dst                 destination pool (initialized by constructor function)
        seed_object         new seed object

        NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
              you should NOT call it when lock can be used by another thread.
        ************************************************************************/
        public static void ae_shared_pool_set_seed(shared_pool dst, alglib.apobject seed_object)
        {
            dst.seed_object = seed_object.make_copy();
            dst.recycled_objects = null;
            dst.enumeration_counter = null;
        }


        /************************************************************************
        This  function  retrieves  a  copy  of  the seed object from the pool and
        stores it to target variable.

        pool                pool
        obj                 target variable
        
        NOTE: this function IS thread-safe.  It  acquires  pool  lock  during its
              operation and can be used simultaneously from several threads.
        ************************************************************************/
        public static void ae_shared_pool_retrieve<T>(shared_pool pool, ref T obj) where T : alglib.apobject
        {
            alglib.apobject new_obj;
            
            /* assert that pool was seeded */
            alglib.ap.assert(pool.seed_object!=null, "ALGLIB: shared pool is not seeded, PoolRetrieve() failed");
            
            /* acquire lock */
            ae_acquire_lock(pool.pool_lock);
            
            /* try to reuse recycled objects */
            if( pool.recycled_objects!=null )
            {
                /* retrieve entry/object from list of recycled objects */
                sharedpoolentry result = pool.recycled_objects;
                pool.recycled_objects = pool.recycled_objects.next_entry;
                new_obj = result.obj;
                result.obj = null;
                
                /* move entry to list of recycled entries */
                result.next_entry = pool.recycled_entries;
                pool.recycled_entries = result;
                
                /* release lock */
                ae_release_lock(pool.pool_lock);
                
                /* assign object to smart pointer */
                obj = (T)new_obj;
                
                return;
            }
                
            /*
             * release lock; we do not need it anymore because
             * copy constructor does not modify source variable.
             */
            ae_release_lock(pool.pool_lock);
            
            /* create new object from seed */
            new_obj = pool.seed_object.make_copy();
                
            /* assign object to pointer and return */
            obj = (T)new_obj;
        }


        /************************************************************************
        This  function  recycles object owned by the source variable by moving it
        to internal storage of the shared pool.

        Source  variable  must  own  the  object,  i.e.  be  the only place where
        reference  to  object  is  stored.  After  call  to  this function source
        variable becomes NULL.

        pool                pool
        obj                 source variable

        NOTE: this function IS thread-safe.  It  acquires  pool  lock  during its
              operation and can be used simultaneously from several threads.
        ************************************************************************/
        public static void ae_shared_pool_recycle<T>(shared_pool pool, ref T obj) where T : alglib.apobject
        {
            sharedpoolentry new_entry;
            
            /* assert that pool was seeded */
            alglib.ap.assert(pool.seed_object!=null, "ALGLIB: shared pool is not seeded, PoolRecycle() failed");
            
            /* assert that pointer non-null */
            alglib.ap.assert(obj!=null, "ALGLIB: obj in ae_shared_pool_recycle() is NULL");
            
            /* acquire lock */
            ae_acquire_lock(pool.pool_lock);
            
            /* acquire shared pool entry (reuse one from recycled_entries or malloc new one) */
            if( pool.recycled_entries!=null )
            {
                /* reuse previously allocated entry */
                new_entry = pool.recycled_entries;
                pool.recycled_entries = new_entry.next_entry;
            }
            else
            {
                /*
                 * Allocate memory for new entry.
                 *
                 * NOTE: we release pool lock during allocation because new() may raise
                 *       exception and we do not want our pool to be left in the locked state.
                 */
                ae_release_lock(pool.pool_lock);
                new_entry = new sharedpoolentry();
                ae_acquire_lock(pool.pool_lock);
            }
            
            /* add object to the list of recycled objects */
            new_entry.obj = obj;
            new_entry.next_entry = pool.recycled_objects;
            pool.recycled_objects = new_entry;
            
            /* release lock object */
            ae_release_lock(pool.pool_lock);
            
            /* release source pointer */
            obj = null;
        }


        /************************************************************************
        This function clears internal list of  recycled  objects,  but  does  not
        change seed object managed by the pool.

        pool                pool

        NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
              you should NOT call it when lock can be used by another thread.
        ************************************************************************/
        public static void ae_shared_pool_clear_recycled(shared_pool pool)
        {
            pool.recycled_objects = null;
        }


        /************************************************************************
        This function allows to enumerate recycled elements of the  shared  pool.
        It stores reference to the first recycled object in the smart pointer.

        IMPORTANT:
        * in case target variable owns non-NULL value, it is rewritten
        * recycled object IS NOT removed from pool
        * target variable DOES NOT become owner of the new value; you can use
          reference to recycled object, but you do not own it.
        * this function IS NOT thread-safe
        * you SHOULD NOT modify shared pool during enumeration (although you  can
          modify state of the objects retrieved from pool)
        * in case there is no recycled objects in the pool, NULL is stored to obj
        * in case pool is not seeded, NULL is stored to obj

        pool                pool
        obj                 reference
        ************************************************************************/
        public static void ae_shared_pool_first_recycled<T>(shared_pool pool, ref T obj) where T : alglib.apobject
        {   
            /* modify internal enumeration counter */
            pool.enumeration_counter = pool.recycled_objects;
            
            /* exit on empty list */
            if( pool.enumeration_counter==null )
            {
                obj = null;
                return;
            }
            
            /* assign object to smart pointer */
            obj = (T)pool.enumeration_counter.obj;
        }


        /************************************************************************
        This function allows to enumerate recycled elements of the  shared  pool.
        It stores pointer to the next recycled object in the smart pointer.

        IMPORTANT:
        * in case target variable owns non-NULL value, it is rewritten
        * recycled object IS NOT removed from pool
        * target pointer DOES NOT become owner of the new value
        * this function IS NOT thread-safe
        * you SHOULD NOT modify shared pool during enumeration (although you  can
          modify state of the objects retrieved from pool)
        * in case there is no recycled objects left in the pool, NULL is stored.
        * in case pool is not seeded, NULL is stored.

        pool                pool
        obj                 target variable
        ************************************************************************/
        public static void ae_shared_pool_next_recycled<T>(shared_pool pool, ref T obj) where T : alglib.apobject
        {   
            /* exit on end of list */
            if( pool.enumeration_counter==null )
            {
                obj = null;
                return;
            }
            
            /* modify internal enumeration counter */
            pool.enumeration_counter = pool.enumeration_counter.next_entry;
            
            /* exit on empty list */
            if( pool.enumeration_counter==null )
            {
                obj = null;
                return;
            }
            
            /* assign object to smart pointer */
            obj = (T)pool.enumeration_counter.obj;
        }


        /************************************************************************
        This function clears internal list of recycled objects and  seed  object.
        However, pool still can be used (after initialization with another seed).

        pool                pool
        state               ALGLIB environment state

        NOTE: this function is NOT thread-safe. It does not acquire pool lock, so
              you should NOT call it when lock can be used by another thread.
        ************************************************************************/
        public static void ae_shared_pool_reset(shared_pool pool)
        {   
            pool.seed_object = null;
            pool.recycled_objects = null;
            pool.enumeration_counter = null;
        }
    }
}
public partial class alglib
{
    /*
     * Parts of alglib.smp class which are NOT shared with GPL version of ALGLIB
     */
    public partial class smp
    {
        #pragma warning disable 420
        
        public const int AE_SMP_MAXPARAMS = 32;
        public const int AE_QUEUE_SIZE = 1024;
        public const int AE_WAIT_CYCLES = 32768;
        public const int AE_QUICK_SCANS_COUNT = 1024;
        public const int AE_SLEEP_ON_IDLE = 1;
        public const int AE_SLEEP_ON_FULL_QUEUE = 1;
        public const int AE_WRK_DISPOSE = 1;
        public const int AE_WRK_NEXT = 2;
        
        /*
         * DESCRIPTION: cores_count  = maximum number of active workers.
         *              cores_to_use = recommended number of active workers:
         *              * positive value >=1 is used to specify exact number of active workers
         *              * 0 means that ALL available cores are used
         *              * negative value means that all cores EXCEPT for cores_to_use will be used
         *                (say, -1 means that all cores except for one will be used.) At least one
         *                core will be used in this case, even if you assign -9999999 to this field.
         * PROTECTION:  not needed
         */
        #if AE_NWORKERS_1
        public static volatile int cores_count = 1;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_2
        public static volatile int cores_count = 2;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_3
        public static volatile int cores_count = 3;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_4
        public static volatile int cores_count = 4;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_5
        public static volatile int cores_count = 5;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_6
        public static volatile int cores_count = 6;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_7
        public static volatile int cores_count = 7;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_8
        public static volatile int cores_count = 8;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_9
        public static volatile int cores_count = 9;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_10
        public static volatile int cores_count = 10;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_11
        public static volatile int cores_count = 11;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_12
        public static volatile int cores_count = 12;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_13
        public static volatile int cores_count = 13;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_14
        public static volatile int cores_count = 14;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_15
        public static volatile int cores_count = 15;
        public static volatile int cores_to_use = 0;
        #elif AE_NWORKERS_16
        public static volatile int cores_count = 16;
        public static volatile int cores_to_use = 0;
        #else
        public static volatile int cores_count = System.Environment.ProcessorCount;
        public static volatile int cores_to_use = -1;
        #endif
        
        public static void AE_CRITICAL_ASSERT(bool x)
        {
            if( !x )
                System.Environment.FailFast("ALGLIB: critical error in thread manager");
        }

        //
        // Delegate which is called from the worker thread.
        //
        public delegate void task_func(ae_task_data data);
        

        /*************************************************************************
        Event.

        This class is a wrapper for NET event.
             
        IMPORTANT: although NET provides fully functional event class, it  is  not
                   recommended to rely on ability of an auto reset event to  hadle
                   multiple threads waiting for it to signal. The reason  is  that
                   C++ version of ae_event class does not support such functionality,
                   and relying on it in the NET version may  lead  to  portability
                   issues.
        *************************************************************************/
        public class ae_event
        {
            public System.Threading.EventWaitHandle event_handle;
        }


        /*************************************************************************
        Thread.

        This class is a wrapper for thread handle
        *************************************************************************/
        public class ae_thread_handle
        {
            public System.Threading.Thread thread;
        }


        /*************************************************************************
        Class which is used to store parameters of the task.
        *************************************************************************/
        public class ae_task_parameter
        {
            public double dval;
            public int ival;
            public bool bval;
            public alglib.complex cval;
            public object val;
        }


        /*************************************************************************
        Task data: set of task parameters and pointer to task function.
        *************************************************************************/
        public class ae_task_data
        {
            public ae_task_parameter[] parameters;
            public task_func func;
        };


        /*************************************************************************
        Task information: task data and synchronization structures
        *************************************************************************/
        public class ae_task_info
        {
            /*
             * DESCRIPTION: auto reset event which is set when root task is solved.
             *              non-signalling by default.
             */
            public ae_event done_event;

            /*
             * DESCRIPTION: task data, undefined by default.
             * PROTECTION:  no lock-based protection. However, it is assumed that
             *              a) before push_task() only one thread works with data
             *              b) push_task() generates memory barrier
             *              c) after pop_task() or steal_task() only one thread
             *                 works with data
             *              These assumptions should give enough of protection.
             */
            public volatile ae_task_data data;

            /*
             * DESCRIPTION: parent task group. can be null. null by default.
             *              when isn't null, task completion is signalled by:
             *              * acquiring parent_group.group_lock
             *              * decreasing waiting_count (and setting exception if needed)
             *              * in case new value of waiting_count is not 0 OR wake_up_worker_on_completion
             *                is false, we release group_lock and exit.
             *              * otherwise, we wake up sleeping worker and
             *                give him our queue by:
             *                * updating worker_idx
             *                * releasing group_lock
             *                * clearing wake_up_worker_on_completion field
             *                * signalling parent_worker.wakeup_event
             *                Finally, we dispose themseles (dispose_thread()
             *                is called and we wait for the wakeup_event)
             *              This field is NULL by default.
             *
             * PROTECTION:  same as for task_info.data
             */
            public volatile ae_task_group parent_group;

            /*
             * DESCRIPTION: a list of child task groups owned by this task.
             *              Groups are linked using ae_task_group->next_group field.
             *              This field is NULL on initial creation.
             *              Every time we create group it is added to this list.
             *              Every time we dispose group, it is removed from the list.
             *
             * PROTECTION:  no protection. Only one worker thread works with this list.
             */
            public volatile ae_task_group child_groups;

            /*
             * DESCRIPTION: set to non-null value when task (or one of its childs)
             *              spawned exception during execution.
             *              null by default.
             * PROTECTION:  no protection. This field is modified in the context
             *              of the thread which works with task.
             */
            public volatile System.Exception exception;

            /*
             * DESCRIPTION: used to organize child tasks into linked list
             *              (or to organize disposed tasks into list).
             * PROTECTION:  not needed (only parent worker references it).
             */
            public volatile ae_task_info next_task;
        }


        /*************************************************************************
        Group of child tasks
        *************************************************************************/
        public class ae_task_group
        {
            /*
             * DESCRIPTION: parent thread which spawned this specific group of tasks.
             *              can NOT be null.
             * PROTECTION:  not needed. This field is modified only by
             *              create_task_group() method which generates appropriate barriers
             *              so its value is visible to all other threads. Other methods
             *              and/or threads can read this value, but can not modify it.
             */
            public volatile ae_worker_thread parent_worker;

            /*
             * DESCRIPTION: number of child tasks waiting for processing (i.e. tasks which
             *              were attached to group and push'ed to queue). This field is
             *              incremented by push_task(), decremented by solve_task().
             * PROTECTION:  by group_lock.
             */
            public volatile int waiting_count;
            
            /*
             * DESCRIPTION: worker thread sets this field to true when it wants to be
             *              waked up on completion of this group. This field is set to
             *              false by thread which wakes up worker.
             * PROTECTION:  by group_lock.
             */
            public volatile bool wake_up_worker_on_completion;

            /*
             * DESCRIPTION: lock which protects group
             * PROTECTION:  by itself
             */
            public ae_lock group_lock;

            /*
             * DESCRIPTION: list of child tasks, updated by push_task().
             *              contains tasks which were not processed immediately
             *              after push_task() in the context of the worker thread.
             * PROTECTION:  not needed (only parent_worker works with list)
             */
            public volatile ae_task_info childs;

            /*
             * DESCRIPTION: set to non-NULL value when child task raised exception 
             *              during execution. NULL by default.
             * PROTECTION:  before group is completed - protected by group_lock.
             *              after completion - it is possible to access this field without locking.
             */
            public volatile System.Exception exception;

            /*
             * DESCRIPTION: this field is used to organize groups in a linked list:
             *              either list of child groups of the worker thread or
             *              list of disposed groups. null by default.
             * PROTECTION:  not needed (only parent worker can modify this list).
             */
            public volatile ae_task_group next_group;
        }


        /*************************************************************************
        Worker thread
        *************************************************************************/
        public class ae_worker_thread
        {
            /*
             * Auto-reset event which is set to signalled state when all child tasks
             * are done (thread, other than the worker thread, which decreases
             * waiting_count down to zero, sets this event).
             *
             * When this event signalled, all tasks were processed (solved
             * or failed), waiting_count is zero, failed_count contains
             * number of failed tasks.
             */
            public ae_event wakeup_event;

            /*
             * DESCRIPTION: thread instance, set during initial creation,
             *              remains unmodified through all the lifetime.
             * PROTECTION:  not needed, other threads do not reference it.
             */
            public ae_thread_handle handle;

            /*
             * DESCRIPTION: worker_idx - index of the worker thread, from 1 to  queue_count-1.
             *              This index is set by create_worker(). It can be modified by another
             *              worker thread which wakes up current worker.
             *
             *              It is responsibility of the party who changed worker_idx to call
             *              ae_pin_thread()/ae_unpin_thread() in order to change information
             *              about ownership. The only exception from this rule is initial
             *              creation, when ownership on the queue is acquired in the worker_loop()
             *              function.
             *
             * PROTECTION:  protected, but has no dedicated lock. Protection protocol
             *              is described below:
             *
             *              1. modification of worker_idx means that thread is assigned
             *                 to some queue or reassigned to new queue
             *
             *              2. worker_idx is modified only when:
             *                 a) thread is created
             *                 b) thread starts to wait for the unfinished child group
             *                     (worker_idx is set to -1)
             *                 c) other thread which finished child group awoke thread
             *                 d) thread is disposed
             *                 e) disposed thread is reused
             *
             *              3. in case (a) worker_idx is set before thread is created,
             *                 protection is provided by implicit barriers which OS
             *                 generates during thread creation
             *
             *              4. in cases (b) and (c) modification of worker_idx is protected
             *                 by group_lock of the group being waited for:
             *                 * in case (b) thread modifies its own worker_idx (sets to -1)
             *                   after acquiring group_lock of group G1 (group it waits for).
             *                 * in case (c) another thread mofifies worker_idx of the thread
             *                   after acquiring group_lock of G1.
             *                 * in both cases worker_idx is modified simultaneously with
             *                   wake_up_worker_on_completion flag
             *                 In these cases protection is provided by the fact that only
             *                 one group can have wake_up_worker_on_completion flag set to
             *                 TRUE, only one thread can finish this group, and all operations
             *                 with this flag and worker_idx are performed within group_lock.
             *
             *              5. in cases (d) and (e) it is guaranteed that all child groups were
             *                 finished prior to disposing/waking up thread. So, no one can
             *                 modify worker_idx in attempt to wakeup thread.
             *                 In these cases protection is provided by thread's wakeup_event:
             *                 * in case (d) worker_idx is set to -1 before waiting for event
             *                 * in case (e) worker_idx is modified before thread wakes up.
             *
             */
            public volatile int worker_idx;

            /*
             * DESCRIPTION: this field is used to organize disposed threads into linked list.
             * PROTECTION:  not needed (only dispose_thread can modify this list).
             */
            public volatile ae_worker_thread next_thread;
        }


        /*************************************************************************
        Task queue
        *************************************************************************/
        public class ae_worker_queue
        {
            /*
             * DESCRIPTION: queue status:
             *              * queue_lock    - lock which protects status fields
             *                                (implemented using interlocked
             *                                operations, 1 when acquired,
             *                                0 when released)
             *              * tasks         - circular buffer of tasks, unused
             *                                elements are equal to null.
             *                                tasks are pushed to top, popped
             *                                from top, stealed from bottom.
             *              * top           - index of top element. Tasks are
             *                                stored in tasks[top], tasks[top+1], ..
             *              * cnt           - number of tasks in a queue
             *              * queue_size    - size of the queue
             * PROTECTION:  by queue_lock.
             */
            public ae_lock queue_lock;
            public volatile ae_task_info[] tasks;
            public volatile int top;
            public volatile int cnt;
            public volatile int queue_size;
        }


        /*************************************************************************
        Thread pool
        *************************************************************************/
        public class ae_thread_pool
        {
            /*
             * DESCRIPTION: queues, including primary queue and worker queues.
             *              can be NULL when queues_count==0.
             * PROTECTION:  not needed (initialized during creation, not changed since then)
             */
            public volatile ae_worker_queue[] queues;
            
            /*
             * DESCRIPTION: total number of queues, including primary queue
             *              and worker queues, >=2. Equal to number of cores+1.
             * PROTECTION:  not needed (initialized during creation, not changed since then)
             */
            public volatile int queues_count;
    
            /*
             * DESCRIPTION: this pair of objects is used to track status of the root tasks.
             *              Every time we push root task, root_cnt is increased. Every time
             *              we solve root task (one with no task group), root_cnt is decreased.
             *
             *              Every time root_cnt becomes nonzero, root_tasks_are_present
             *              event is set to signalling state. Every time root_cnt becomes
             *              zero, root_tasks_are_present event is set to non-signalling.
             *
             *              ae_push_root_task() is responsible for increase of root_cnt and
             *              setting event to signalling, ae_solve_task() is responsible for
             *              decrease of root_cnt and clearing event.
             *
             * PROTECTION:  both fields are protected by queues[0].queue_lock.  No protection
             *              when queues_count==0. Although events have their own protection,
             *              we MUST set/unset event only when lock is acquired.
             */
            public volatile int root_cnt;
            public ae_event root_tasks_are_present;
            
            /*
             * DESCRIPTION: whether termination request was submitted or not
             * PROTECTION:  by queues[0].queue_lock
             *              no protection when queues_count==0.
             */
            public volatile bool termination_required;
            
            /*
             * DESCRIPTION: pool of disposed tasks
             * PROTECTION:  tasks_lock
             */
            public volatile ae_task_info disposed_tasks;
            public ae_lock tasks_lock;
            
            /*
             * DESCRIPTION: pool of disposed groups
             * PROTECTION:  groups_lock
             */
            public volatile ae_task_group disposed_groups;
            public ae_lock groups_lock;
            
            /*
             * DESCRIPTION: pool of disposed worker threads
             * PROTECTION:  threads_lock
             */
            public volatile ae_worker_thread disposed_threads;
            public ae_lock threads_lock;
        }
        
    
        /*
         * ALGLIB thread pool
         */
        public static ae_thread_pool main_thread_pool;

        //
        // Task-specific information:
        // * tsk_current_worker - contains null when we work from the
        //   external thread, non-null when task is solved from the
        //   worker thread.
        // * tsk_current_task - null by default, initialized by ae_solve_task()
        //   when we start working on task, restored to previous value on exit
        //   from ae_solve_task. Similar actions are performed by ae_push_task()
        //   and ae_push_root_task() when we solve tasks immediately in the
        //   context of the current thread.
        //
        [System.ThreadStatic]
        public static ae_worker_thread  tsk_current_worker;
        
        [System.ThreadStatic]
        public static ae_task_info      tsk_current_task;

        /*
         * Static constructor.
         *
         * NOTE: it is translation of the original C function, so it
         *       contains a lot of unneeded initializations by nulls or zeros.
         *       However, it is better not to change code which works, so
         *       all these unneeded initializations are copied as is.
         */
        static smp()
        {
            int i, j;
            
            /*
             * Determine cores count.
             * Do not initialize pool on single-core systems.
             */
            if( cores_count<2 )
                return;
    
            /*
             * Allocate
             */
            main_thread_pool = new ae_thread_pool();
    
            /*
             * Initialize pool object as if we have no multicore support.
             */
            main_thread_pool.queues_count = 0;
            main_thread_pool.queues = null;
            main_thread_pool.root_cnt = 0;
            main_thread_pool.termination_required = false;
            main_thread_pool.disposed_threads = null;
            main_thread_pool.disposed_tasks = null;
            main_thread_pool.disposed_groups = null;
            ae_init_event(ref main_thread_pool.root_tasks_are_present, true);
            ae_init_lock(ref main_thread_pool.threads_lock);
            ae_init_lock(ref main_thread_pool.tasks_lock);
            ae_init_lock(ref main_thread_pool.groups_lock);
    
            /*
             * Initialize queues
             */
            main_thread_pool.queues_count = cores_count + 1;
            main_thread_pool.queues = new ae_worker_queue[main_thread_pool.queues_count];
            for(i=0; i<main_thread_pool.queues_count; i++)
            {
                main_thread_pool.queues[i] = new ae_worker_queue();
                ae_init_lock(ref main_thread_pool.queues[i].queue_lock);
                main_thread_pool.queues[i].top = 0;
                main_thread_pool.queues[i].cnt = 0;
                main_thread_pool.queues[i].queue_size = AE_QUEUE_SIZE;
                main_thread_pool.queues[i].tasks = new ae_task_info[AE_QUEUE_SIZE];
                for(j=0; j<AE_QUEUE_SIZE; j++)
                    main_thread_pool.queues[i].tasks[j] = null;
            }
            
            
            /*
             * Create worker threads
             */
            for(i=1; i<main_thread_pool.queues_count; i++)
                ae_create_worker(i);
        }


        /************************************************************************
        This function returns number of CPU cores which should be used by  worker
        threads, as specified by user. In case user specified non-positive number
        of  cores  to  use,  this numner will be converted according to following
        rules:
        *  0 => cores_count
        * -1 => max(cores_count-1,1)
        * -2 => max(cores_count-2,1)
        and so on.

        This function requires initialized thread pool. It will  fail  if  called
        without initialized thread pool.
        ************************************************************************/
        public static int ae_cores_to_use()
        {
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            if( cores_to_use<=0 )
            {
                int r = (main_thread_pool.queues_count-1)+cores_to_use;
                return r>=1 ? r : 1;
            }
            else
                return cores_to_use;
        }

        /************************************************************************
        This function initializes ae_event structure and sets to non-signalling
        mode.
        ************************************************************************/
        public static void ae_init_event(ref ae_event obj, bool manual_reset)
        {
            obj = new ae_event();
            if( manual_reset )
                obj.event_handle = new System.Threading.ManualResetEvent(false);
            else
                obj.event_handle = new System.Threading.AutoResetEvent(false);
        }


        /************************************************************************
        This function waits for event
        ************************************************************************/
        public static void ae_wait_for_event(ae_event obj)
        {
            obj.event_handle.WaitOne();
        }


        /************************************************************************
        This function sets event to signalling state
        ************************************************************************/
        public static void ae_set_event(ae_event obj)
        {
            obj.event_handle.Set();
        }


        /************************************************************************
        This function sets event to nonsignalling state
        ************************************************************************/
        public static void ae_reset_event(ae_event obj)
        {
            obj.event_handle.Reset();
        }


        /************************************************************************
        This function frees ae_event structure.
        ************************************************************************/
        public static void ae_free_event(ref ae_event obj)
        {
            obj = null;
        }

        /************************************************************************
        This function starts the thread and stores its handle  into  thread_info.
        It provides OS-independent abstraction layer for thread creation.

        PARAMETERS:
            thread_function     -   main function
            instance            -   instance of ae_worker_thread, must be fully
                                    inialized by the caller except for
                                    instance->handle parameter which is initialized
                                    by this function

        NOTE: after  its  creation new thread MUST wait for instace->wakeup_event
              which will be set to signalling state after thread handle  will  be
              successfully stored in the thread instance structure.
              
        NOTE: this function should NOT be called when AE_OS is AE_UNKNOWN  -  the
              whole program will be abnormally terminated.
        ************************************************************************/
        public static void ae_start_thread(
            System.Threading.ParameterizedThreadStart thread_function,
            ae_worker_thread instance)
        {
            instance.handle.thread = new System.Threading.Thread(thread_function);
            instance.handle.thread.IsBackground = true;
            instance.handle.thread.Priority = System.Threading.ThreadPriority.Highest;
            ae_set_event(instance.wakeup_event);
            instance.handle.thread.Start(instance);
        }

        /************************************************************************
        This function exits from the current thread.
        It provides OS-independent abstraction layer for thread finalization.
              
        NOTE: this function should NOT be called when AE_OS is AE_UNKNOWN  -  the
              whole program will be abnormally terminated.
        ************************************************************************/
        public static void ae_exit_thread(ae_worker_thread instance)
        {
            instance.handle.thread.Abort();
        }


        /************************************************************************
        This function pauses current thread for specified number of milliseconds.
        It provides OS-independent abstraction layer for Sleep() call.
              
        NOTE: this function should NOT be called when AE_OS is AE_UNKNOWN  -  the
              whole program will be abnormally terminated.
        ************************************************************************/
        public static void ae_sleep_thread(int ms_to_sleep)
        {
            System.Threading.Thread.Sleep(ms_to_sleep);
        }


        /************************************************************************
        This function sets worker_idx of the thread and (depending on  OS)  tries
        to pin thread to its personal core.
        
        NOTE: thread->worker_idx must be clear (negative) prior to  calling  this
              function.
              
        NOTE: no synchronization is used during this call, it is your responsibility
              to ensure that access to worker_idx is synchronized.
        ************************************************************************/
        public static void ae_pin_thread(ae_worker_thread thread, int worker_idx)
        {
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(thread.worker_idx<0);
            AE_CRITICAL_ASSERT(worker_idx>0 && worker_idx<main_thread_pool.queues_count);
            thread.worker_idx = worker_idx;
            //SetThreadAffinityMask(thread->thread_handle.thread, 1L<<(worker_idx-1));
        }
        
        
        /************************************************************************
        This function clears thread->worker_idx and clears thread  affinity  mask
        of the worker.
        ************************************************************************/
        public static void ae_unpin_thread(ae_worker_thread thread)
        {
            AE_CRITICAL_ASSERT(thread.worker_idx>0);
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            thread.worker_idx = -1;
            //SetThreadAffinityMask(thread->thread_handle.thread, (1L<<(main_thread_pool->queues_count-1))-1);
        }


        /************************************************************************
        This function  creates  worker  thread  and  assigns it to specific queue
        given by worker_idx. New worker thread is active and running  immediately
        after its creation.

        It is assumed that this method is called  by  the  worker   thread  which
        previously worked with that queue (owns the queue). After this  call, new
        worker thread owns the queue, and previous worker has lost its ownership.

        NOTE: worker_idx>=1
              
        NOTE: this function should NOT be called when AE_OS is AE_UNKNOWN  -  the
              whole program will be abnormally terminated.
        ************************************************************************/
        public static void ae_create_worker(int worker_idx)
        {
            ae_worker_thread thread;
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(worker_idx>0 && worker_idx<main_thread_pool.queues_count);

            /* thread object: use list of disposed threads or start new thread */
            ae_acquire_lock(main_thread_pool.threads_lock);
            if( main_thread_pool.disposed_threads==null )
            {
                /* release lock */
                ae_release_lock(main_thread_pool.threads_lock);

                /* no disposed threads, create new one */
                thread = new ae_worker_thread();
                ae_init_event(ref thread.wakeup_event, false);
                thread.handle = new ae_thread_handle();
                thread.worker_idx = worker_idx;
                thread.next_thread = null;
            
                /* start thread; new thread is automatically pinned to appropriate core in the worker_loop function. */ 
                ae_start_thread(ae_worker_loop, thread);
            }
            else
            {
                /* we have thread in the pool, use it */
                thread = main_thread_pool.disposed_threads;
                main_thread_pool.disposed_threads = thread.next_thread;
                
                /* release lock */
                ae_release_lock(main_thread_pool.threads_lock);
                
                /* initialize fields */
                AE_CRITICAL_ASSERT(thread.worker_idx<0);
                thread.next_thread = null;
                
                /* pin thread to appropriate core (it is VERY important to pin thread BEFORE waking it up) */
                ae_pin_thread(thread, worker_idx);
                
                /* wake up thread */
                ae_set_event(thread.wakeup_event);
            }
        }


        /************************************************************************
        This function disposes worker thread -  worker thread is pauses and moved
        to internal list of paused workers which can be reused later.

        NOTE: it is expected that thread to be  disposed  has  no  child  groups.
              critical error is raised when list has unfinished childs.
              
        NOTE: this function should NOT be called when we have no threading support
              - the whole program will be terminated.
        ************************************************************************/
        public static void ae_dispose_worker(ae_worker_thread thread)
        {
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(thread.worker_idx<0); 
            
            /* move thread to list of disposed threads, wait for wakeup */
            ae_acquire_lock(main_thread_pool.threads_lock);
            thread.next_thread = main_thread_pool.disposed_threads;
            main_thread_pool.disposed_threads = thread;
            ae_release_lock(main_thread_pool.threads_lock);
            ae_wait_for_event(thread.wakeup_event);
        }
        
        
        /*************************************************************************
        This function solves task, handles errors.

        For non-root tasks it reports results to task group, for root ones  -   it
        decreases root_cnt field and sets event to signalling state.

        Boolean parameter from_queue is true when task was extracted  from  queue,
        false when task was executed without pushing it to queue.  The  difference
        is that in the second case we do NOT decrement parent_group->waiting_count.

        Return value:
        * AE_WRK_NEXT in case task was successfully solved and we have  to  search
          for the next one.
        * AE_WRK_DISPOSE in case task was successfully solved and waked  up  other
          worker thread which now owns our queue. Current worker  thread  have  to
          dispose itself.
          
        The only situation when AE_WRK_DISPOSE can be returned is  when  there was
        some worker thread (different from current worker) which  posted  task  to
        the queue and called wait_for_group(). Thus, if:
        * from_queue is false => AE_WRK_NEXT is returned
        * task was posted by current worker => AE_WRK_NEXT is returned

        NOTE 1: _state parameter stores information about threading state (current
                worker thread and current thread pool).
                
        NOTE 2: It is expected that this function is called only from  the  worker
                thread. Any attempt to call it from the program main  thread  will
                lead to abnormal termination of the program.
                
        NOTE 3: this function solves task and makes sure that no unprocessed child
                groups left. It terminates program in case  unprocessed  group  is
                detected.
                
        NOTE:   on exception this function sets task->exception AND
                task->parent_group->exception to current instance of exception.
                Exception is NOT rethrown - only silently signalled.
        *************************************************************************/
        public static int ae_solve_task(ae_task_info task, bool from_queue)
        {
            ae_worker_thread sleeping_worker;
            ae_task_info prev_task;
            System.Exception exception;
            
            AE_CRITICAL_ASSERT(task!=null);
            AE_CRITICAL_ASSERT(task.child_groups==null);
            
            /*
             * Create state structure, solve task
             */
            prev_task = tsk_current_task;
            tsk_current_task = task;
            try
            {
                exception = null;
                try
                {
                    task.data.func(task.data);
                }
                catch(System.Exception e)
                {
                    exception = e;
                    ae_terminate_child_groups();
                }
                AE_CRITICAL_ASSERT(task.child_groups==null);
                task.exception = exception;
            }
            finally
            {
                tsk_current_task = prev_task;
            }
            
            /*
             * Problem is solved, postprocessing
             */
            if( task.parent_group!=null )
            {
                /*
                 * task is a part of some group:
                 * 1. decrease waiting count for this group
                 * 2. update group->exception
                 * 2. in case some thread waits for completion, wake it up and give our queue to this thread
                 */
                ae_task_group parent_group = task.parent_group;
                if( from_queue )
                {
                    ae_acquire_lock(parent_group.group_lock);
                    parent_group.waiting_count--;
                    if( exception!=null )
                        parent_group.exception = exception;
                    if( parent_group.waiting_count==0 && parent_group.wake_up_worker_on_completion )
                    {
                        int queue_idx;
                        
                        /*
                         * There is some worker thread which waits for completion of this group.
                         * We wake up this thread and give our queue to it.
                         *
                         * NOTE: this branch of code is NOT executed when we work without SMP support.
                         *       In this case tasks are executed immediately when they are pushed to
                         *       queue - BEFORE parent task calls wait_for_group() and sets 
                         *       wake_up_worker_on_completion field.
                         *       However, we perform several safety checks here.
                         */
                        AE_CRITICAL_ASSERT(main_thread_pool!=null);
                        AE_CRITICAL_ASSERT(main_thread_pool.queues_count>=2);
                        AE_CRITICAL_ASSERT(tsk_current_worker!=null);
                        AE_CRITICAL_ASSERT(parent_group.parent_worker!=null);
                        AE_CRITICAL_ASSERT(parent_group.parent_worker!=tsk_current_worker);
                        AE_CRITICAL_ASSERT(parent_group.parent_worker.worker_idx<0);  
                        
                        /* cache variables */
                        queue_idx = tsk_current_worker.worker_idx;
                        sleeping_worker = parent_group.parent_worker;
                        
                        /* change group state, unpin one worker, pin another one */ 
                        parent_group.wake_up_worker_on_completion = false;
                        ae_unpin_thread(tsk_current_worker); 
                        ae_pin_thread(sleeping_worker, queue_idx); 
                        
                        /* release lock, notify sleeping worker that it is time to wake up */ 
                        ae_release_lock(parent_group.group_lock);
                        ae_set_event(sleeping_worker.wakeup_event);
                        
                        return AE_WRK_DISPOSE;
                    }
                    else
                        ae_release_lock(parent_group.group_lock);
                }
                else
                {
                    if( exception!=null )
                    {
                        ae_acquire_lock(parent_group.group_lock);
                        parent_group.exception = exception;
                        ae_release_lock(parent_group.group_lock);
                    }
                }
                return AE_WRK_NEXT;
            }
            else
            {
                /*
                 * We've solved root task.
                 *
                 * NOTE: this branch of code is NOT executed when we work without SMP support.
                 *       In this case root tasks are solved in push_root_task() without calling
                 *       this function. However, we perform several safety checks here.
                 */
                AE_CRITICAL_ASSERT(main_thread_pool!=null);
                ae_acquire_lock(main_thread_pool.queues[0].queue_lock);
                main_thread_pool.root_cnt--;
                if( main_thread_pool.root_cnt==0 )
                    ae_reset_event(main_thread_pool.root_tasks_are_present);
                ae_release_lock(main_thread_pool.queues[0].queue_lock);
                ae_set_event(task.done_event);
                return AE_WRK_NEXT;
            }
        }


        /*************************************************************************
        This  function  creates  new  instance of task_info structure, either root
        task (in case no parent group is specified)  or  child  task  attached  to
        group passed as argument.

        Additional parameter _state must point to the current state of the  ALGLIB
        environment. Among other information it stores current threading settings.
        Checks listed in notes 2-3  are  performed  using  information  stored  in
        _state.

        NOTE 1: all fields of the task_info object have their default  values,  as
                specified in the task_info description.
                
        NOTE 2: in case parent_group is NULL, we must call this function from  the
                external (non-worker) thread. Only non-worker  thread  can  create
                root task.
                
        NOTE 3: in case parent_group is non-NULL, we must call this function  from
                the worker thread. External threads  can  not  create  groups  and
                attach tasks to them.
        *************************************************************************/
        public static ae_task_info ae_create_task(ae_task_group parent_group)
        {
            ae_task_info task;
            int i;
            
            /* quick exit for system without SMP support */
            if( main_thread_pool==null )
            {
                task = new ae_task_info();
                ae_init_event(ref task.done_event, false);
                task.exception = null;
                task.parent_group = parent_group;
                task.child_groups = null;
                if( parent_group!=null )
                {
                    task.next_task = parent_group.childs;
                    parent_group.childs = task;
                }
                else
                    task.next_task = null;
                task.data = new ae_task_data();
                task.data.parameters = new ae_task_parameter[AE_SMP_MAXPARAMS];
                for(i=0; i<AE_SMP_MAXPARAMS; i++)
                    task.data.parameters[i] = new ae_task_parameter();
                return task;
            }
            
            /* allocate memory, primary initialization */
            ae_acquire_lock(main_thread_pool.tasks_lock);
            if( main_thread_pool.disposed_tasks==null )
            {
                /* release lock */
                ae_release_lock(main_thread_pool.tasks_lock);

                /* no disposed tasks, create new one */
                task = new ae_task_info();
                ae_init_event(ref task.done_event, false);
                task.data = new ae_task_data();
                task.data.parameters = new ae_task_parameter[AE_SMP_MAXPARAMS];
                for(i=0; i<AE_SMP_MAXPARAMS; i++)
                    task.data.parameters[i] = new ae_task_parameter();
            }
            else
            {
                /* we have thread in the pool, use it */
                task = main_thread_pool.disposed_tasks;
                main_thread_pool.disposed_tasks = task.next_task;
                
                /* release lock */
                ae_release_lock(main_thread_pool.tasks_lock);
            }
            
            /* initialization of other fields */
            task.exception = null;
            task.parent_group = parent_group;
            task.child_groups = null;
            if( parent_group!=null )
            {
                task.next_task = parent_group.childs;
                parent_group.childs = task;
            }
            else
                task.next_task = null;

            /* exit */
            return task;
        }


        /*************************************************************************
        This function disposes instance of ae_task_info structure by  freeing  all
        dynamically allocated structures. After call to this function reference to
        ae_task_info becomes invalid.

        This  function  may store structure in the internal list for reuse.
        *************************************************************************/
        public static void ae_dispose_task(ae_task_info task)
        {
            /* check correctness of parameters */
            AE_CRITICAL_ASSERT(task!=null);
            AE_CRITICAL_ASSERT(task.child_groups==null);
            
            /* dispose task depending on SMP support */
            if( main_thread_pool==null )
            {
                /*
                 * quick exit for OS without SMP support
                 */
                ae_free_event(ref task.done_event);
                task = null;
            }
            else
            {
                /*
                 * OS support for SMP is detected.
                 * Move task to list of disposed tasks
                 */
                ae_acquire_lock(main_thread_pool.tasks_lock);
                task.next_task = main_thread_pool.disposed_tasks;
                main_thread_pool.disposed_tasks = task;
                ae_release_lock(main_thread_pool.tasks_lock);
            }
        }


        /*************************************************************************
        This function creates new  instance  of  task_group.  The  task  group  is
        attached to the parent task (as specified by _state->parent_task).

        Additional parameter _state must point to the current state of the  ALGLIB
        environment. Among other information it stores current threading settings.

        NOTE: this function may be called from non-worker thread, in this case  it
              returns NULL. It also returns NULL when OS provides no  support  for
              SMP (main_thread_pool==NULL).
        *************************************************************************/
        public static ae_task_group ae_create_task_group()
        {
            ae_task_group group;
            
            /*
             * No SMP support is present
             */
            if( main_thread_pool==null )
            {
                group = new ae_task_group();
                ae_init_lock(ref group.group_lock);
                AE_CRITICAL_ASSERT(tsk_current_task!=null);
                group.parent_worker = null;
                group.waiting_count = 0;
                group.wake_up_worker_on_completion = false;
                group.childs = null;
                group.exception = null;
                group.next_group = tsk_current_task.child_groups;
                tsk_current_task.child_groups = group;
                return group;
            }
            
            /* allocate memory, primary initialization */
            ae_acquire_lock(main_thread_pool.groups_lock);
            if( main_thread_pool.disposed_groups==null )
            {
                /* release lock */
                ae_release_lock(main_thread_pool.groups_lock);

                /* no disposed groups, create new one */
                group = new ae_task_group();
                ae_init_lock(ref group.group_lock);
            }
            else
            {
                /* we have thread in the pool, use it */
                group = main_thread_pool.disposed_groups;
                main_thread_pool.disposed_groups = group.next_group;
                
                /* release lock */
                ae_release_lock(main_thread_pool.groups_lock);
            }
            
            /* initialize other fields */
            AE_CRITICAL_ASSERT(tsk_current_task!=null);
            group.parent_worker = tsk_current_worker;
            group.waiting_count = 0;
            group.wake_up_worker_on_completion = false;
            group.childs = null;
            group.exception = null;
            group.next_group = tsk_current_task.child_groups;
            tsk_current_task.child_groups = group;

            /* exit */
            return group;
        }


        /*************************************************************************
        This function disposes instance of ae_task_group structure by  freeing all
        dynamically allocated structures. After call to this function reference to
        ae_task_group becomes invalid.

        This function may store structure in the internal list for reuse.

        NOTE: all  tasks  must  be  disposed  prior  to  calling this function. It
              simply ignores presence of child tasks.
              
        NOTE: this function can be used with NULL group.
        *************************************************************************/
        public static void ae_dispose_group(ae_task_group group)
        {
            /* dispose group depending on SMP support */
            if( main_thread_pool==null )
            {
                /*
                 * free memory - for OS without SMP support
                 */
                ae_free_lock(ref group.group_lock);
                group = null;
            }
            else
            {
                /*
                 * Move group to list of disposed groups
                 */
                if( group==null )
                    return;
                ae_acquire_lock(main_thread_pool.groups_lock);
                group.next_group = main_thread_pool.disposed_groups;
                main_thread_pool.disposed_groups = group;
                ae_release_lock(main_thread_pool.groups_lock);
            }
        }


        /*************************************************************************
        This  function  waits for the completion of the task group  owned  by  the
        current task. Only task which owns the group can call this method.

        This method:
        * tries to execute child tasks in the context of  the  current  worker  by
          traversing list of the child tasks,  calling ae_pop_specific_task()  and
          executing them. It is important to traverse list of the  child  problems
          in the correct direction - starting from the most recently added task.
        * if all child tasks were executed in the context of  the current  worker,
          we return.
        * if one of the child tasks was not found on top of the  stack,  it  means
          that it was stolen. In this case we wait for the group to  complete  and
          give our queue to other worker thread.  After  waking  up  (possibly  in
          another queue) we return.
          
        NOTE 1: this  function  can  be  called  from  non-worker thread with NULL
                group.

        NOTE 2: when called with group=null, this function silently  returns  back
                to  caller.  Group  can  be  NULL when called from any kind of the
                thread.

        NOTE 3: depending on dispose_on_success parameter function may either:
                * dispose all child tasks AND group itself (group is removed  from
                  childs of the current task)
                * dispose all child tasks, but leave group in the default  (empty)
                  state. Group is not removed from the childs of the task.
                
        NOTE 4: this function rethrows exceptions raised in child tasks. In case of
                multiple exceptions only one is rethrown.

        NOTE 5: on exception this function disposes only group and its child tasks,
                but leaves other child groups of the current task unchanged. It is
                responsibility of ae_solve_task to dispose other groups.

        NOTE:   this method may reassign current worker object to another queue
                in case we had to actually wait for the childs to complete.
        *************************************************************************/
        public static void ae_wait_for_group(
            ae_task_group group,
            bool dispose_on_success)
        {
            ae_task_info task;
            System.Exception exception;
            int queue_idx;
            
            /* check consistency of the threading information */
            AE_CRITICAL_ASSERT( group==null || group.parent_worker==tsk_current_worker );
            
            /* quick exit */
            if( group==null )
                return;
            
            /* start waiting for childs */
            ae_acquire_lock(group.group_lock);
            if( group.waiting_count>0 )
            {
                /* 
                 * There are childs waiting for processing.
                 * Try to process them within context of the current thread.
                 *
                 * NOTE: this branch of code is executed only when we have SMP support because
                 *       in other cases tasks are executed immediately when they are push'ed.
                 *       However, we perform several safety checks here.
                 *
                 * NOTE: we do not protect access to group->childs because only owner of the
                 *       group may modify its childs - and we are the owner
                 */
                AE_CRITICAL_ASSERT( main_thread_pool!=null );
                AE_CRITICAL_ASSERT( tsk_current_worker!=null );
                AE_CRITICAL_ASSERT( group!=null );
                AE_CRITICAL_ASSERT( group.childs!=null );
                ae_release_lock(group.group_lock);
                for(task = group.childs;
                    (task!=null) && (ae_pop_specific_task(task,tsk_current_worker.worker_idx)!=null);
                    task = task.next_task)
                {
                    int tmp;
                    
                    tmp = ae_solve_task(task, false);
                    AE_CRITICAL_ASSERT(tmp==AE_WRK_NEXT);
                    ae_acquire_lock(group.group_lock);
                    AE_CRITICAL_ASSERT(group.waiting_count>0);
                    group.waiting_count--;
                    ae_release_lock(group.group_lock);
                }
                
                /* in case there are still exist unprocessed childs,
                   wait for them to be completed by other threads */
                ae_acquire_lock(group.group_lock);
                if( group.waiting_count>0 )
                {
                    group.wake_up_worker_on_completion = true;
                    queue_idx = tsk_current_worker.worker_idx;
                    ae_unpin_thread(tsk_current_worker);
                    ae_release_lock(group.group_lock);
                    ae_create_worker(queue_idx);
                    ae_wait_for_event(tsk_current_worker.wakeup_event);
                    AE_CRITICAL_ASSERT(tsk_current_worker.worker_idx>0);
                    ae_acquire_lock(group.group_lock);
                    exception = group.exception;
                    ae_release_lock(group.group_lock);
                }
                else
                {
                    exception = group.exception;
                    ae_release_lock(group.group_lock);
                }
            }
            else
            {
                exception = group.exception;
                ae_release_lock(group.group_lock);
            }
                
            /* childs are done, dispose child tasks */
            while( group.childs!=null )
            {
                task = group.childs;
                group.childs = task.next_task;
                ae_dispose_task(task);
            }
            
            /* dispose group itself (if needed) */
            if( exception!=null || dispose_on_success )
            {
                /* remove group from list of child groups */
                if( tsk_current_task.child_groups!=group )
                {
                    /* group to be removed is not first in the list */
                    for(ae_task_group cg = tsk_current_task.child_groups; cg!=null; cg=cg.next_group)
                        if( cg.next_group==group )
                        {
                            cg.next_group = group.next_group;
                            break;
                        }
                }
                else
                    tsk_current_task.child_groups = group.next_group;
                                
                /* group is removed from list, we can dispose it */
                ae_dispose_group(group);
            }
            
            /* rethrow exception if needed */
            if( exception!=null )
                throw exception;
        }


        /*************************************************************************
        This  function  terminates child groups of parent task  (as  specified  by
        the tsk_current_task structure).

        This function must be called prior to raising exception in order  to  make
        sure  that  all  child  tasks  are  done  before  we  start  to deallocate
        dynamically allocated data structures.

        NOTE: this function does not terminate parent task itself.
        *************************************************************************/
        public static void ae_terminate_child_groups()
        {
            AE_CRITICAL_ASSERT( tsk_current_task!=null );
            
            while( tsk_current_task.child_groups!=null )
                ae_terminate_group(tsk_current_task.child_groups);
        }

        /*************************************************************************
        This  function  terminates  task  group  by  removing its tasks from queue
        (tasks are not solved - just removed) and waiting for completion of  tasks
        which were stolen by other workers. Only task which  owns  the  group  can
        call this method. Finally, we dispose group.

        NOTE:   this function does not rethrow exceptions in child tasks.
        *************************************************************************/
        public static void ae_terminate_group(ae_task_group group)
        {
            ae_task_info task;
            int queue_idx;
            
            /* check consistency of the threading information */
            AE_CRITICAL_ASSERT( tsk_current_task!=null );
            AE_CRITICAL_ASSERT( group==null || group.parent_worker==tsk_current_worker );
            
            /* quick exit */
            if( group==null )
                return;
            
            /* start waiting for childs */
            ae_acquire_lock(group.group_lock);
            if( group.waiting_count>0 )
            {
                /* 
                 * There are childs waiting for processing.
                 * Try to remove them without solving.
                 *
                 * NOTE: this branch of code is executed only when we have SMP support because
                 *       in other cases tasks are executed immediately when they are push'ed.
                 *       However, we perform several safety checks here.
                 *
                 * NOTE: we do not protect access to group->childs because only owner of the
                 *       group may modify its childs - and we are the owner
                 */
                AE_CRITICAL_ASSERT( main_thread_pool!=null );
                AE_CRITICAL_ASSERT( tsk_current_worker!=null );
                AE_CRITICAL_ASSERT( group.childs!=null );
                ae_release_lock(group.group_lock);
                for(task = group.childs;
                    (task!=null) && (ae_pop_specific_task(task,tsk_current_worker.worker_idx)!=null);
                    task = task.next_task)
                {
                    ae_acquire_lock(group.group_lock);
                    AE_CRITICAL_ASSERT(group.waiting_count>0);
                    group.waiting_count--;
                    ae_release_lock(group.group_lock);
                }
                
                /* in case there are still exist unprocessed childs,
                   wait for them to be completed by other threads */
                ae_acquire_lock(group.group_lock);
                if( group.waiting_count>0 )
                {
                    group.wake_up_worker_on_completion = true;
                    queue_idx = tsk_current_worker.worker_idx;
                    ae_unpin_thread(tsk_current_worker); 
                    ae_release_lock(group.group_lock);
                    ae_create_worker(queue_idx);
                    ae_wait_for_event(tsk_current_worker.wakeup_event);
                    AE_CRITICAL_ASSERT(tsk_current_worker.worker_idx>0);
                }
                else
                    ae_release_lock(group.group_lock);
            }
            else
                ae_release_lock(group.group_lock);
                
            /* childs are done, dispose child tasks */
            while( group.childs!=null )
            {
                task = group.childs;
                group.childs = task.next_task;
                ae_dispose_task(task);
            }
            
            /* dispose group itself */
            if( tsk_current_task.child_groups!=group )
            {
                /* group to be removed is not first in the list */
                for(ae_task_group cg = tsk_current_task.child_groups; cg!=null; cg=cg.next_group)
                    if( cg.next_group==group )
                    {
                        cg.next_group = group.next_group;
                        break;
                    }
            }
            else
                tsk_current_task.child_groups = group.next_group;
            ae_dispose_group(group);
        }


        /************************************************************************
        This  function  pushes  non-root  task  object  to the top of the current
        worker queue.

        NOTE 1: this method must act as full barrier, i.e. all pending
                modifications to task_info instance will be committed to the
                memory before return from this method.

        NOTE 2: in case queue has free space, task is added to the queue. In case
                queue is full, this task is solved immediately without modyfing queue.

        NOTE 3: _state parameter stores information about threading state (current
                worker thread and current thread pool).
                
        NOTE 4: It is expected that this function is called only from the child
                task. Any attempt to call it from the program main thread will
                lead to abnormal termination of the program.

        NOTE 5: task MUST be part of some task group. Any attempt to use this
                function on a task which is not part of the group will terminate
                the program.
                
        NOTE 6: parent group MUST be owned by the same worker as one which calls
                push_task(). Attempts to push tasks tied to groups owned by other
                workers will crash program.
        ************************************************************************/
        public static void ae_push_task(ae_task_info task, bool execute_immediately)
        {
            ae_worker_queue queue;
            
                
            /*
             * Handle situation when no threading support is present or
             * task must be executed immediately
             */
            if( main_thread_pool==null || execute_immediately )
            {
                int tmp;
                AE_CRITICAL_ASSERT(tsk_current_task!=null);
                tmp = ae_solve_task(task, false);
                AE_CRITICAL_ASSERT(tmp==AE_WRK_NEXT);
                return;
            }
            
            /*
             * Threading support is present. Check parameters.
             */
            AE_CRITICAL_ASSERT(tsk_current_task!=null);
            AE_CRITICAL_ASSERT(tsk_current_worker!=null);
            AE_CRITICAL_ASSERT(task!=null);
            AE_CRITICAL_ASSERT(task.parent_group!=null);
            AE_CRITICAL_ASSERT(task.parent_group.parent_worker==tsk_current_worker);
            
            /*
             * Try to push task to queue
             */
            queue = main_thread_pool.queues[tsk_current_worker.worker_idx];
            ae_acquire_lock(queue.queue_lock);
            if( queue.cnt>=queue.queue_size )
            {
                int tmp;
                ae_release_lock(queue.queue_lock);
                tmp = ae_solve_task(task, false);
                AE_CRITICAL_ASSERT(tmp==AE_WRK_NEXT);
            }
            else
            {
                queue.top = queue.top==0 ? queue.queue_size-1 : queue.top-1;
                queue.cnt++;
                queue.tasks[queue.top] = task;
                ae_release_lock(queue.queue_lock);
                ae_acquire_lock(task.parent_group.group_lock);
                task.parent_group.waiting_count++;
                ae_release_lock(task.parent_group.group_lock);
            }
        }


        /*************************************************************************
        This function pushes root task object to  the  top   of   the  main  queue
        (one with index 0). After task added to queue this  function  returns.  If
        you  want  to  wait  for  the  task  completion,  you  have  to  wait  for
        task->done_event.

        NOTE 1: this method must act as full barrier, i.e. all pending
                modifications to task_info instance will be committed to the
                memory before return from this method.

        NOTE 2: in case queue has free space, task is added to the queue. In  case
                queue is full, we wait for some time before trying add something.
                
        NOTE 3: It is expected that this function is called  only  from  the  main
                thread. Any attempt to call it from the worker  thread  will  lead
                to abnormal termination of the program.

        NOTE 4: task MUST must NOT be part of some task group. Any attempt to  use
                this function on a task which is part of the group will  terminate
                the program.
                
        NOTE 5: this method correctly handles situation when  we  have  NO  worker
                threads by solving task in the context of the calling thread. Such
                situation is possible, for example, when AE_OS==AE_UNKNOWN
        *************************************************************************/
        public static void ae_push_root_task(ae_task_info task)
        {
            ae_worker_queue queue;
            AE_CRITICAL_ASSERT(task!=null);
            AE_CRITICAL_ASSERT(task.parent_group==null);
            
            /*
             * Handle situation when no SMP support is detected
             */
            if( main_thread_pool==null )
            {
                ae_task_info prev_task = tsk_current_task;
                tsk_current_task = task;
                try
                {
                    System.Exception exception = null;
                    try
                    {
                        task.data.func(task.data);
                    }
                    catch(System.Exception e)
                    {
                        exception = e;
                        ae_terminate_child_groups();
                    }
                    AE_CRITICAL_ASSERT(task.child_groups==null);
                    task.exception = exception;
                }
                finally
                {
                    tsk_current_task = prev_task;
                }
                ae_set_event(task.done_event);
                return;
            }
            
            /*
             * SMP support is present, post task to queue
             */
            AE_CRITICAL_ASSERT(main_thread_pool.queues_count>=2);
            queue = main_thread_pool.queues[0];
            ae_acquire_lock(queue.queue_lock);
            for(;;)
            {
                /*
                 * pause execution in case queue is full
                 */
                if( queue.cnt==queue.queue_size)
                {
                    ae_release_lock(queue.queue_lock);
                    ae_sleep_thread(AE_SLEEP_ON_FULL_QUEUE); /* TODO: better way of pausing execution */
                    ae_acquire_lock(queue.queue_lock);
                    continue;
                }
                
                /*
                 * we have free space, insert task and return
                 */
                queue.tasks[(queue.top+queue.cnt)%queue.queue_size] = task;
                queue.cnt++;
                if( main_thread_pool.root_cnt==0 )
                    ae_set_event(main_thread_pool.root_tasks_are_present);
                main_thread_pool.root_cnt++;
                ae_release_lock(queue.queue_lock);
                return;
            }        
        }


        /*************************************************************************
        This method pops task from the top of the  corresponding  queue.  In  case
        queue has items, this method will return non-NULL pointer. In  case  queue
        is empty, it will return NULL.

        NOTE 1: this method can pop elements from any queue, independently of  the
                queue owner - oven from queues which belong to  other  threads  or
                from the root queue.
        *************************************************************************/
        public static ae_task_info ae_pop_task(int queue_idx)
        {
            ae_worker_queue queue;
            
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(queue_idx>=0 && queue_idx<main_thread_pool.queues_count);
            
            queue = main_thread_pool.queues[queue_idx];
            ae_acquire_lock(queue.queue_lock);
            if( queue.cnt==0 )
            {
                ae_release_lock(queue.queue_lock);
                return null;
            }
            else
            {
                ae_task_info task = queue.tasks[queue.top];
                queue.tasks[queue.top] = null;
                queue.cnt--;
                queue.top = (queue.top+1)%queue.queue_size;
                ae_release_lock(queue.queue_lock);
                return task;
            }
        }


        /*************************************************************************
        This method tries to pop specific task from the top of the   corresponding
        queue. In case specific task is not found in the queue, it returns NULL.

        NOTE 1: this method can pop elements from any queue, independently of  the
                queue owner - oven from queues which belong to  other  threads  or
                from the root queue.
        *************************************************************************/
        public static ae_task_info ae_pop_specific_task(ae_task_info task, int queue_idx)
        {
            ae_worker_queue queue;
            
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(queue_idx>=0 && queue_idx<main_thread_pool.queues_count);
            
            queue = main_thread_pool.queues[queue_idx];
            ae_acquire_lock(queue.queue_lock);
            if( queue.cnt==0 || queue.tasks[queue.top]!=task )
            {
                ae_release_lock(queue.queue_lock);
                return null;
            }
            else
            {
                queue.tasks[queue.top] = null;
                queue.cnt--;
                queue.top = (queue.top+1)%queue.queue_size;
                ae_release_lock(queue.queue_lock);
                return task;
            }
        }

        /*************************************************************************
        This method steals task from the BOTTOM of the  corresponding  queue.   In
        case queue has items, this method will return non-NULL pointer.  In   case
        queue is empty, it will return NULL.

        NOTE 1: this method can steal elements from any  queue,  independently  of
                the queue owner - oven from queues which belong to  other  threads
                or from the root queue.
        *************************************************************************/
        public static ae_task_info ae_steal_task(int queue_idx)
        {
            ae_worker_queue queue;
            
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(queue_idx>=0 && queue_idx<main_thread_pool.queues_count);
            
            queue = main_thread_pool.queues[queue_idx];
            ae_acquire_lock(queue.queue_lock);
            if( queue.cnt==0 )
            {
                ae_release_lock(queue.queue_lock);
                return null;
            }
            else
            {
                int idx;
                ae_task_info task;
                idx = (queue.top+queue.cnt-1)%queue.queue_size;
                task = queue.tasks[idx];
                queue.tasks[idx] = null;
                queue.cnt--;
                ae_release_lock(queue.queue_lock);

                /* exit */
                return task;
            }
        }


        /*************************************************************************
        This method is used to update SMP status of ALGLIB function.

        Every SMP-capable ALGLIB function can work in multi-threaded  and  single-
        threaded  modes.  When  working  in  single-threaded  mode, SPAWN operator
        executes  task  immediately  in  the context of the current worker thread.
        When  working  in  multithreaded  mode,  SPAWN operator creates child task
        which is attached to automatically created group and pushed to queue.

        ALGLIB  functions  can switch between two modes during their execution. In
        order  to  store  information  about  current  mode, each function has two
        automatically created variables:
        * _child_tasks  -   group of child tasks, which is NULL  by  default,  but
                            automatically created when SMP support is turned on.
        * _smp_enabled  -   current SMP status, false by default.

        When SMP is enabled, task_group object is automatically created.  It  will
        never  be  freed  until  the  end  of  the  current  function. When SMP is
        disabled,  child  group  is  left non-NULL. Thus,  three  combinations  of
        variable values are possible:
        *  _smp_enabled && _child_tasks!=NULL       SMP is active, group created
        * !_smp_enabled && _child_tasks!=NULL       SMP is inactive, but was previously
                                                    active; we can have child tasks
                                                    in the queue.
        * !_smp_enabled && _child_tasks==NULL       SMP was inactive since the beginning
                                                    of the current function.

        This method  updates  SMP  status  and  changes  state  of  the  variables
        according to new state and their  previous  values.  It  is  important  to
        understand, that this method changes only values of  variables  passed  by
        reference.  It  does  NOT  changes  some  global  settings, like number of
        working threads, threads status and so on.

        NOTE: this function does not change anything  when  called from non-worker
              thread. SMP support can be activated only  from  the  worker  thread
              executed as part of the ALGLIB thread pool.
              When it is called from external, non-worker thread, it just silently
              returns without changing SMP status.
        *************************************************************************/
        public static void ae_set_smp_support(
            ref ae_task_group _child_tasks,
            ref bool _smp_enabled,
            bool new_support_status)
        {
            if( main_thread_pool==null || tsk_current_worker==null )
            {
                AE_CRITICAL_ASSERT(_child_tasks==null && !_smp_enabled);
                return;
            }
            if( new_support_status )
            {
                if( !_smp_enabled )
                {
                    if( _child_tasks==null )
                        _child_tasks = ae_create_task_group();
                    _smp_enabled = true;
                }
                AE_CRITICAL_ASSERT(_child_tasks!=null);
            }
            else
                _smp_enabled = false;
        }

        /*************************************************************************
        This method is used  to  synchronize  with  child  problems  according  to
        presence of child tasks and current status of SMP support.

        Every SMP-capable ALGLIB function can work in multi-threaded  and  single-
        threaded modes. It is possible to switch between two modes and to be in  a
        single-threaded  mode,  but  to have child tasks previously spawned in the
        multi-threaded mode.

        This function performs synchronization with childs:
        * in case code is executed in context of the  non-worker  thread,  it just
          silently returns. We check that no task group was created and SMP is off
          and abort program in case these conditions  are  not  satisfied  (it  is
          considered to be critical error in the program logic).
        * in case we have NULL  task  group,  we  silently  return.  Similarly  to
          previous case we check that SMP support is turned off.
          
        *************************************************************************/
        public static void ae_sync(
            ae_task_group _child_tasks,
            bool smp_enabled,
            bool dispose_group)
        {
            if( main_thread_pool==null || tsk_current_worker==null )
            {
                AE_CRITICAL_ASSERT(_child_tasks==null && !smp_enabled);
                return;
            }
            if( _child_tasks==null )
            {
                AE_CRITICAL_ASSERT(!smp_enabled);
                return;
            }
            ae_wait_for_group(_child_tasks, dispose_group);
        }
        
        
        public static void ae_worker_loop(System.Object obj)
        {
            ae_worker_thread thread;
            int result, hint_index, queue_idx;
            
            thread = (ae_worker_thread)obj;
            
            /*
             * Wait for wakeup_event
             * Pin thread to worker queue (some tricky ops with worker_idx are
             * required because ae_pin_thread needs thread with clear worked_idx).
             */
            ae_wait_for_event(thread.wakeup_event);
            queue_idx = thread.worker_idx;
            thread.worker_idx = -1;
            ae_pin_thread(thread, queue_idx);
    
            /*
             * Outer loop:
             * 1. execute inner loop (until exit from loop)
             * 2. wait for pool->root_tasks_are_present
             * 3. goto (1)
             */
            tsk_current_worker = thread;
            AE_CRITICAL_ASSERT(thread!=null);
            AE_CRITICAL_ASSERT(main_thread_pool!=null);
            AE_CRITICAL_ASSERT(thread.worker_idx>0 && thread.worker_idx<main_thread_pool.queues_count);
            for(;;)
            {
                /*
                 * Inner loop:
                 * 0. in case worker_idx is greater than number of cores to use (as requested by user),
                 *    sleep for 1 ms, then exit inner loop
                 * 1. scan for new tasks to process: pop from our queue, steal
                 *    from other queues. In case task is found, goto (4).
                 *    Goto (2) otherwise.
                 * 2. in case pool->root_cnt is zero, exit inner loop
                 * 3. wait for AE_WAIT_CYCLES cycles, goto (0)
                 * 4. solve task
                 * 5. depending on task result, dispose worker (we awoke another
                 *    worker which owns our queue) or goto (0)
                 */
                hint_index = -1;
                for(;;)
                {
                    ae_task_info current_task;
                    int i, j;
                    AE_CRITICAL_ASSERT(thread.worker_idx>0 && thread.worker_idx<main_thread_pool.queues_count);
            
                    /*
                     * (0) check that worker index is less than or equal to number of cores which are allowed to use.
                     *     This block allows user to dynamically control number of cores which can be utilized
                     *     by ALGLIB.
                     *
                     *     In case this specific worker is restricted from performing activity, we perform short sleep
                     *     for AE_SLEEP_ON_IDLE ms, and then we jump back to outer loop.
                     *
                     *     NOTE 1: short period of sleep is important because it prevents us from utilizing all CPU
                     *             power while continuously checking for permission to work.
                     *
                     *     NOTE 2: worker may leave some unprocessed tasks in its queue; it is not problem because
                     *             these tasks will be stolen by some other worker.
                     */
                    if( thread.worker_idx>ae_cores_to_use() )
                    {
                        ae_sleep_thread(AE_SLEEP_ON_IDLE);
                        break;
                    }
                    
                    /*
                     * (1) scan for new tasks
                     *     * first, we try to pop from our own queue
                     *     * second, we try to steal task from quueue hint_index (in case hint_index is non-negative).
                     *       in any case hint_index is reset to -1 after this block.
                     *     * then, we try to steal from non-root queues
                     *     * finally, we try to steal from root queue  
                     */
                    AE_CRITICAL_ASSERT(hint_index<main_thread_pool.queues_count);
                    current_task = ae_pop_task(thread.worker_idx);
                    if( current_task==null && hint_index>=0 )
                        current_task = ae_steal_task(hint_index);
                    hint_index = -1;
                    if( current_task==null )
                    {
                        int offs = thread.worker_idx;
                        int cnt = main_thread_pool.queues_count;
                        for(i=0; i<=cnt; i++)
                            if( (i+offs)%cnt!=thread.worker_idx && (i+offs)%cnt!=0 )
                            {
                                current_task = ae_steal_task((i+offs)%cnt);
                                if( current_task!=null )
                                    break;
                            }
                    }
                    if( current_task==null )
                        current_task = ae_steal_task(0);
                    
                    /*
                     * (2), (3): no task found
                     */
                    if( current_task==null )
                    {
                        /*
                         * No tasks found.
                         * Depending on presense of active root tasks we either 
                         * a) BREAK to outer cycle (no root tasks running)
                         * b) wait and CONTINUE to beginning of inner cycle (active inner task is running)
                         *    hint_index variable is set to index of non-empty queue in case
                         *    such queue is found during waiting.
                         */
                        bool no_root_tasks;
                        
                        /* breaking to ourer cycle if needed */
                        ae_acquire_lock(main_thread_pool.queues[0].queue_lock);
                        no_root_tasks = main_thread_pool.root_cnt==0;
                        ae_release_lock(main_thread_pool.queues[0].queue_lock);
                        if( no_root_tasks )
                            break;
                        
                        /* wait (we scan queues during waiting) and continue inner cycle */
                        hint_index = -1;
                        for(i=0; i<=AE_QUICK_SCANS_COUNT; i++)
                        {
                            for(j=0; j<main_thread_pool.queues_count; j++)
                                if( main_thread_pool.queues[j].cnt>0 )
                                    hint_index = j;
                            if( hint_index>=0 )
                                break;
                        }
                        if( hint_index<0 )
                            ae_yield();
                        continue;
                    }
                    
                    /*
                     * (4) execute task
                     */
                    AE_CRITICAL_ASSERT(thread.worker_idx>=0 && thread.worker_idx<main_thread_pool.queues_count);
                    result = ae_solve_task(current_task, true);
                    AE_CRITICAL_ASSERT( (result==AE_WRK_DISPOSE) || (thread.worker_idx>=0 && thread.worker_idx<main_thread_pool.queues_count) );
                    
                    // TODO: check that after successful solution of the task worker thread have no child groups.
                    //       in case there exists unprocessed group, user-mode exception should be generated.
                    //       the only situation where there can be unprocessed groups is when task generated exception
                    //       during its execution.
                    
                    /*
                     * (5) postprocessing
                     */
                    if( result==AE_WRK_DISPOSE )
                        ae_dispose_worker(thread);
                }
                
                /*
                 * We've exited from inner loop, it means that no active root
                 * tasks is present. Wait for new tasks to come.
                 */
                ae_wait_for_event(main_thread_pool.root_tasks_are_present);
            }
        }
    }
    
    /*************************************************************************
    SMP self-tests.
    *************************************************************************/
    public class smpselftests
    {
        /*************************************************************************
        Test 1: test for basic processing and exception handling.

        Task parameters:
        - params[0].val      array to sort
        - params[1].val      temporary buffer
        - params[2].ival     left bound of [i0,i1) subinterval to sort
        - params[3].ival     right bound of [i0,i1) subinterval to sort
        - params[4].ival     task with i0=params[4].ival=i1-1 throws exception
        *************************************************************************/

        /* test 1: child function - merge sorted subarrays */
        public static void test1_merge_func(smp.ae_task_data data)
        {
            int idx0, idx1, idx2, srcleft, srcright, dst;
            int[] a, buf;

            a    = (int[])data.parameters[0].val;
            buf  = (int[])data.parameters[1].val;
            idx0 = data.parameters[2].ival;
            idx1 = data.parameters[3].ival;
            idx2 = data.parameters[4].ival;
            
            srcleft = idx0;
            srcright = idx1;
            dst = idx0;
            for(;;)
            {
                if( srcleft==idx1&&srcright==idx2 )
                {
                    break;
                }
                if( srcleft==idx1 )
                {
                    buf[dst] = a[srcright];
                    srcright = srcright+1;
                    dst = dst+1;
                    continue;
                }
                if( srcright==idx2 )
                {
                    buf[dst] = a[srcleft];
                    srcleft = srcleft+1;
                    dst = dst+1;
                    continue;
                }
                if( a[srcleft]<a[srcright] )
                {
                    buf[dst] = a[srcleft];
                    srcleft = srcleft+1;
                    dst = dst+1;
                }
                else
                {
                    buf[dst] = a[srcright];
                    srcright = srcright+1;
                    dst = dst+1;
                }
            }
            for(dst=idx0; dst<=idx2-1; dst++)
                a[dst] = buf[dst];
        }

        /* root function for test 1 */
        public static void test1_root_func(smp.ae_task_data data)
        {
            /*
             * unload parameters
             */
            int i0, i1, idxa, eidx;
            int[] arr, buf;
            smp.ae_task_group group0, group1;
            smp.ae_task_info  task0, task1;
            
            arr  = (int[])data.parameters[0].val;
            buf  = (int[])data.parameters[1].val;
            i0   = data.parameters[2].ival;
            i1   = data.parameters[3].ival;
            eidx = data.parameters[4].ival;
            
            /* exit on unit subproblem */
            if( i0==i1-1 )
            {
                alglib.ap.assert(eidx<i0 || eidx>=i1, "exception generated");
                return;
            }
            
            /* split subproblem into [i0, idxa) and [idxa, i1) */
            do { idxa = i0+math.randominteger(i1-i0); } while( idxa==i0 || idxa==i1 );
            group0 = smp.ae_create_task_group();
            task0 = smp.ae_create_task(group0);
            task0.data.func = test1_root_func;
            task0.data.parameters[0].val = arr;
            task0.data.parameters[1].val = buf;
            task0.data.parameters[2].ival = i0;
            task0.data.parameters[3].ival = idxa;
            task0.data.parameters[4].ival = eidx;
            smp.ae_push_task(task0, false);
            task1 = smp.ae_create_task(group0);
            task1.data.func = test1_root_func;
            task1.data.parameters[0].val = arr;
            task1.data.parameters[1].val = buf;
            task1.data.parameters[2].ival = idxa;
            task1.data.parameters[3].ival = i1;
            task1.data.parameters[4].ival = eidx;
            smp.ae_push_task(task1, false);
            smp.ae_wait_for_group(group0, true);
            
            /* merge */
            group1 = smp.ae_create_task_group();
            task0 = smp.ae_create_task(group1);
            task0.data.func = test1_merge_func;
            task0.data.parameters[0].val = arr;
            task0.data.parameters[1].val = buf;
            task0.data.parameters[2].ival = i0;
            task0.data.parameters[3].ival = idxa;
            task0.data.parameters[4].ival = i1;
            smp.ae_push_task(task0, false);
            smp.ae_wait_for_group(group1, true);
        }

        /* this function returns true on success, false on failure */
        public static bool perform_test1()
        {
            int n, pass, i, j, k, eidx;
            int[] arr, buf;
            smp.ae_task_info task;
            bool result;

            n = 100000;
            arr = new int[n];
            buf = new int[n];
            result = true;
            task = null;
            for(pass=0; pass<2; pass++)
            {
                /* decide whether we want to raise exception or to solve problem successfully */
                if( pass==0 )
                    eidx = -1;
                else
                    eidx = math.randominteger(n);
            
                /* create task */
                for(i=0; i<n; i++)
                    arr[i] = i;
                for(i=0; i<n; i++)
                {
                    j = math.randominteger(n);
                    if( j!=i )
                    {
                        k = arr[i];
                        arr[i] = arr[j];
                        arr[j] = k;
                    }
                }
                task = smp.ae_create_task(null);
                task.data.func = test1_root_func;
                task.data.parameters[0].val = arr;
                task.data.parameters[1].val = buf;
                task.data.parameters[2].ival = 0;
                task.data.parameters[3].ival = n;
                task.data.parameters[4].ival = eidx;
                
                /* push task, no exception should be generated at this moment */
                smp.ae_push_root_task(task);
                
                /* wait for the solution, no exception should be generated at this moment */
                smp.ae_wait_for_event(task.done_event);
            
                /* check solution status */
                if( eidx>=0 && eidx<n && task.exception==null )
                    return false;
                if( (eidx<0 || eidx>n) && task.exception!=null )
                    return false;
                if( eidx<0 || eidx>n )
                    for(i=0; i<n; i++)
                        if( arr[i]!=i )
                            return false;
                
                /* dispose */
                smp.ae_dispose_task(task);
                task = null;
            }
            return result;
        }


        /*************************************************************************
        SMP self-tests.

        Test 2: this test spawns A LOT OF root tasks. The idea is to
                check scheduler's ability to handle such stream of tasks.

        Task input parameters:
        - params[0].ival     input value

        Task output parameters:
        - params[0].ival     input*input+1

        *************************************************************************/

        /* root function for test 2 */
        public static void test2_root_func(smp.ae_task_data data)
        {
            int v;
            v = data.parameters[0].ival;
            data.parameters[0].ival = v*v+1;
        }

        /* this function returns true on success, false on failure */
        public static bool perform_test2()
        {
            smp.ae_task_info[] tasks;
            bool result;
            int n, i;

            n = 100000;
            result = true;
            
            /* allocate space */
            tasks = new smp.ae_task_info[n];
            
            /*
             * Create and push tasks.
             * NOTE: we store tasks into the array before pushing in order
             *       to push them as fast as possible.
             */
            for(i=0; i<n; i++)
            {
                tasks[i] = smp.ae_create_task(null);
                tasks[i].data.parameters[0].ival = i;
                tasks[i].data.func = test2_root_func;
            }
            for(i=0; i<n; i++)
                smp.ae_push_root_task(tasks[i]);
            
            /*
             * wait for tasks and check results
             */
            for(i=0; i<n; i++)
            {
                smp.ae_wait_for_event(tasks[i].done_event);
                result = result && (tasks[i].data.parameters[0].ival==i*i+1);
            }
            
            /* dispose */
            for(i=0; i<n; i++)
                smp.ae_dispose_task(tasks[i]);
            
            return result;
        }


        /*************************************************************************
        Test 3: this test spawns A LOT OF child tasks. The idea is to
                check scheduler's ability to handle such stream of tasks.

        Child task input parameters:
        - params[0].ival     input value
        - params[1].val      double* to store result=sin(input)

        Root task input parameters:
        - params[0].ival     N

        Root task output parameters:
        - params[0].dval     sum of child outputs for i=0..N-1
        *************************************************************************/

        /* child function for test 3 */
        public static void test3_child_func(smp.ae_task_data data)
        {
            double v;
            
            /* slow down problem a bit - problems should require
               some amount of time in order to let queue overflow */
            smp.ae_spin_wait(10000);
            
            /* solve problem */
            v = data.parameters[0].ival;
            ((double[])data.parameters[1].val)[data.parameters[2].ival] = System.Math.Sin(v);
        }


        /* root function for test 3 */
        public static void test3_root_func(smp.ae_task_data data)
        {
            smp.ae_task_info[]  tasks;
            smp.ae_task_group[] groups;
            double[] results;
            int n, i;
            double result;

            n = data.parameters[0].ival;
            
            /* allocate space */
            tasks   = new smp.ae_task_info[n];
            groups  = new smp.ae_task_group[n];
            results = new double[n];
            
            /*
             * Create and push tasks.
             * NOTE: we store tasks into the array before pushing in order
             *       to push them as fast as possible.
             */
            for(i=0; i<n; i++)
            {
                groups[i] = smp.ae_create_task_group();
                tasks[i] = smp.ae_create_task(groups[i]);
                tasks[i].data.parameters[0].ival = i;
                tasks[i].data.parameters[1].val  = results;
                tasks[i].data.parameters[2].ival = i;
                tasks[i].data.func = test3_child_func;
                results[i] = 0.0;
            }
            for(i=0; i<n; i++)
                smp.ae_push_task(tasks[i], false);
            
            /*
             * wait for tasks and check results.
             *
             * NOTE: it is important to wait for groups in backward order,
             *       i.e. most recently added group will be waited first.
             *
             *       In fact, it is possible to wait for groups in arbitrary
             *       order, but when we wait for the most recent group, is is
             *       removed from the beginning of the list - almost immediately.
             *       And when we wait for a group in another end of the list,
             *       we have to spent too much time traversing list and removing
             *       group from it.
             */
            result = 0;
            for(i=n-1; i>=0; i--)
            {
                smp.ae_wait_for_group(groups[i], true);
                result += results[i];
            }
            data.parameters[0].dval = result;
        }

        /* this function returns true on success, false on failure */
        public static bool perform_test3()
        {
            smp.ae_task_info task;
            bool result;
            int n, i;
            double t;

            n = 99991;
            task = smp.ae_create_task(null);
            task.data.parameters[0].ival = n;
            task.data.func = test3_root_func;
            smp.ae_push_root_task(task);
            smp.ae_wait_for_event(task.done_event);
            t = 0;
            for(i=0; i<n; i++)
                t = t+System.Math.Sin((double)i);
            result = System.Math.Abs(task.data.parameters[0].dval-t)<=1.0E-6*System.Math.Abs(t);
            smp.ae_dispose_task(task);
            return result;
        }
        

        /*************************************************************************
        SMP self-tests.

        Returns true on success, false on failures.
        *************************************************************************/
        public static bool runtests()
        {
            bool t1, t2, t3;
            
            t1 = perform_test1();
            t2 = perform_test2();
            t3 = perform_test3();
            
            return t1 && t2 && t3;
        }
    }
    
    /************************************************************************
    This function sets number of CPU cores which should  be  used  by  worker
    threads. In case user specified non-positive number of cores to use, this
    number will be converted according to following rules:
    *  0 => ae_cores_count()
    * -1 => max(ae_cores_count()-1,1)
    * -2 => max(ae_cores_count()-2,1)
    and so on.

    In case user specified positive number of  cores,  greater  than 1,  then
    ALGLIB will launch no more than ncores threads (or less, when nworkers is
    larger than actual number of cores).
    ************************************************************************/
    public static void setnworkers(int nworkers)
    {
        alglib.smp.cores_to_use = nworkers;
    }
}
