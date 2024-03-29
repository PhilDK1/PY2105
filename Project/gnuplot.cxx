

/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.c
  @author   N. Devillard
  @date Sep 1998
  @version  $Revision: 2.10 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/
/*--------------------------------------------------------------------------*/

/*
    $Id: gnuplot_i.c,v 2.10 2003/01/27 08:58:04 ndevilla Exp $
    $Author: ndevilla $
    $Date: 2003/01/27 08:58:04 $
    $Revision: 2.10 $
 */

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @file     gnuplot_i.h
  @author   N. Devillard
  @date     Sep 1998
  @version  $Revision: 1.11 $
  @brief    C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.

*/
/*--------------------------------------------------------------------------*/

/*
    $Id: gnuplot_i.h,v 1.11 2003/01/27 08:58:04 ndevilla Exp $
    $Author: ndevilla $
    $Date: 2003/01/27 08:58:04 $
    $Revision: 1.11 $
 */

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

/*---------------------------------------------------------------------------
                                Includes
 ---------------------------------------------------------------------------*/
#include <stdio.h>

/** Maximal number of simultaneous temporary files */
#define GP_MAX_TMP_FILES    64

/*---------------------------------------------------------------------------
                                New Types
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @typedef  gnuplot_ctrl
  @brief    gnuplot session handle (opaque type).

  This structure holds all necessary information to talk to a gnuplot
  session. It is built and returned by gnuplot_init() and later used
  by all functions in this module to communicate with the session, then
  meant to be closed by gnuplot_close().

  This structure is meant to remain opaque, you normally do not need
  to know what is contained in there.
 */
/*-------------------------------------------------------------------------*/

typedef struct _GNUPLOT_CTRL_ {
    /** Pipe to gnuplot process */
    FILE    * gnucmd ;

    /** Number of currently active plots */
    int       nplots ;
    /** Current plotting style */
    char      pstyle[32] ;

    /** Pointer to table of names of temporary files */
    char*      tmp_filename_tbl[GP_MAX_TMP_FILES] ;
    /** Number of temporary files */
    int       ntmp ;
} gnuplot_ctrl ;

/*---------------------------------------------------------------------------
                        Function ANSI C prototypes
 ---------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
/*--------------------------------------------------------------------------*/
gnuplot_ctrl * gnuplot_init(void);

/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
/*--------------------------------------------------------------------------*/
void gnuplot_close(gnuplot_ctrl * handle);

/*-------------------------------------------------------------------------*/
/**
  @brief    Sends a command to an active gnuplot session.
  @param    handle Gnuplot session control handle
  @param    cmd    Command to send, same as a printf statement.

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  @code
  gnuplot_cmd(g, "plot %d*x", 23.0);
  gnuplot_cmd(g, "plot %g * cos(%g * x)", 32.0, -3.0);
  @endcode

  Since the communication to the gnuplot process is run through
  a standard Unix pipe, it is only unidirectional. This means that
  it is not possible for this interface to query an error status
  back from gnuplot.
 */
/*--------------------------------------------------------------------------*/
void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...);

/*-------------------------------------------------------------------------*/
/**
  @brief    Change the plotting style of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_style Plotting-style to use (character string)
  @return   void

  The provided plotting style is a character string. It must be one of
  the following:

  - lines
  - points
  - linespoints
  - impulses
  - dots
  - steps
  - errorbars
  - boxes
  - boxeserrorbars
 */
/*--------------------------------------------------------------------------*/
void gnuplot_setstyle(gnuplot_ctrl * h, char const * plot_style);

/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/
void gnuplot_set_xlabel(gnuplot_ctrl * h, char const * label);


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/
void gnuplot_set_ylabel(gnuplot_ctrl * h, char const * label);

/*-------------------------------------------------------------------------*/
/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
/*--------------------------------------------------------------------------*/
void gnuplot_resetplot(gnuplot_ctrl * h);

/*-------------------------------------------------------------------------*/
/**
  @brief    Plots a 2d graph from a list of doubles.
  @param    handle  Gnuplot session control handle.
  @param    d       Array of doubles.
  @param    n       Number of values in the passed array.
  @param    title   Title of the plot.
  @return   void

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  @code
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot_x(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/
void gnuplot_plot_x(gnuplot_ctrl * handle, double * d, int n, char * title);

/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a 2d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

  @code
    gnuplot_ctrl    *h ;
    double          x[50] ;
    double          y[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        x[i] = (double)(i)/10.0 ;
        y[i] = x[i] * x[i] ;
    }
    gnuplot_plot_xy(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/
void gnuplot_plot_xy(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    char            *   title
) ;


/*-------------------------------------------------------------------------*/
/**
  @brief    Open a new session, plot a signal, close the session.
  @param    title   Plot title
  @param    style   Plot style
  @param    label_x Label for X
  @param    label_y Label for Y
  @param    x       Array of X coordinates
  @param    y       Array of Y coordinates (can be NULL)
  @param    n       Number of values in x and y.
  @return

  This function opens a new gnuplot session, plots the provided
  signal as an X or XY signal depending on a provided y, waits for
  a carriage return on stdin and closes the session.

  It is Ok to provide an empty title, empty style, or empty labels for
  X and Y. Defaults are provided in this case.
 */
/*--------------------------------------------------------------------------*/
/*void gnuplot_plot_once(
    char    *   title,
    char    *   style,
    char    *   label_x,
    char    *   label_y,
    double  *   x,
    double  *   y,
    int         n
);*/

/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a slope on a gnuplot session.
  @param    handle      Gnuplot session control handle.
  @param    a           Slope.
  @param    b           Intercept.
  @param    title       Title of the plot.
  @return   void

  Plot a slope on a gnuplot session. The provided slope has an
  equation of the form y=ax+b

  Example:

  @code
    gnuplot_ctrl    *   h ;
    double              a, b ;

    h = gnuplot_init() ;
    gnuplot_plot_slope(h, 1.0, 0.0, "unity slope") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/
void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
) ;

/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a curve of given equation y=f(x).
  @param    h           Gnuplot session control handle.
  @param    equation    Equation to plot.
  @param    title       Title of the plot.
  @return   void

  Plots out a curve of given equation. The general form of the
  equation is y=f(x), you only provide the f(x) side of the equation.

  Example:

  @code
        gnuplot_ctrl    *h ;
        char            eq[80] ;

        h = gnuplot_init() ;
        strcpy(eq, "sin(x) * cos(2*x)") ;
        gnuplot_plot_equation(h, eq, "sine wave", normal) ;
        gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/
void gnuplot_plot_equation(gnuplot_ctrl * h, char * equation, char * title) ;

/**
 * Writes a CSV file for use with gnuplot commands later.  Allows files to also be saved for post
 * analysis with excel for example. Arguments are similar to gnuplot_plot_x()
 *
 * Uses title as gnuplot "comment" on the first line.
 *
 * @author Peter H Maresh 11/8/2011
 *
 * @param fileName file name to write to.  This should be the same name used in the gnuplot command
 * @param d
 * @param n
 * @param title
 * @return int      <0 if error writing file.
 *
 */
int gnuplot_write_x_csv(
    char const * fileName,
    double const * d,
    int n,
    char const * title);

/**
 * Writes a CSV file for use with gnuplot commands later.  Allows files to also be saved for post
 * analysis with excel for example. Arguments are similar to gnuplot_plot_xy()
 *
 * Uses title as gnuplot "comment" on the first line.
 *
 * @author Peter H Maresh 11/8/2011
 *
 * @param fileName file name to write to.  This should be the same name used in the gnuplot command
 * @param x
 * @param y
 * @param n
 * @param title
 * @return int <0 if file wasn't written.
 */
int gnuplot_write_xy_csv(
    char const *        fileName,
    double const    *   x,
    double const    *   y,
    int                 n,
    char const      *   title);

/**
 * Writes a multi column CSV file for use with gnuplot commands later.  Allows files to also be
 * saved for post analysis with excel for example. Note that when used with gnuplot, since there
 * may be more than one column the whole "1:3" or whatever should be used.
 *
 * Uses title as gnuplot "comment" on the first line.
 *
 * @author Peter H Maresh 11/8/2011
 *
 * @param fileName      file name to write to.  This should be the same name used in the gnuplot
 *                      command
 * @param xListPtr      A list of pointers to column buffers
 * @param n
 * @param numColumns    Length of xListPtr list.  Note that the CSV file will add an additional
 *                      "index" column first.
 * @param title         Title to write for the first line of the .csv file, will be preceeded by
 *                      "#"
 * @return int <0 if file wasn't written.
 */
int gnuplot_write_multi_csv(
    char const *        fileName,
    double const    **  xListPtr,
    int                 n,
    int                 numColumns,
    char const      *   title);

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>

#ifdef WIN32
#include <io.h>
#endif // #ifdef WIN32

/*---------------------------------------------------------------------------
                                Defines
 ---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------
                          Prototype Functions
 ---------------------------------------------------------------------------*/

/**
 * Creates a temporary file name for writing
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 *
 * @return char const * Pointer to file name string.
 */
char const * gnuplot_tmpfile(gnuplot_ctrl * handle);

/**
 * Plot a temporary file.
 *
 * @author Peter (12/9/2011)
 *
 * @param handle
 * @param tmp_filename
 * @param title
 */
void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title);

/*---------------------------------------------------------------------------
                            Function codes
 ---------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
/**
  @brief    Opens up a gnuplot session, ready to receive commands.
  @return   Newly allocated gnuplot control structure.

  This opens up a new gnuplot session, ready for input. The struct
  controlling a gnuplot session should remain opaque and only be
  accessed through the provided functions.

  The session must be closed using gnuplot_close().
 */
/*--------------------------------------------------------------------------*/

gnuplot_ctrl * gnuplot_init(void)
{
    gnuplot_ctrl *  handle ;
    int i;

#ifndef WIN32
    if (getenv("DISPLAY") == NULL) {
        fprintf(stderr, "cannot find DISPLAY variable: is it set?\n") ;
    }
#endif // #ifndef WIN32


    /*
     * Structure initialization:
     */
    handle = (gnuplot_ctrl*)malloc(sizeof(gnuplot_ctrl)) ;
    handle->nplots = 0 ;
    gnuplot_setstyle(handle, (char *) "points") ;
    handle->ntmp = 0 ;
#ifdef WIN32
    char name [200] = "\"C:\\Program Files (x86)\\gnuplot\\bin\\gnuplot.exe\"";
#else
    char name [200] = "gnuplot";
#endif
    
    handle->gnucmd = popen(name, "w") ;
    if (handle->gnucmd == NULL) {
        fprintf(stderr, "error starting gnuplot, is gnuplot or gnuplot.exe in your path?\n") ;
        free(handle) ;
        return NULL ;
    }

    for (i=0;i<GP_MAX_TMP_FILES; i++)
    {
        handle->tmp_filename_tbl[i] = NULL;
    }
    return handle;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Closes a gnuplot session previously opened by gnuplot_init()
  @param    handle Gnuplot session control handle.
  @return   void

  Kills the child PID and deletes all opened temporary files.
  It is mandatory to call this function to close the handle, otherwise
  temporary files are not cleaned and child process might survive.

 */
/*--------------------------------------------------------------------------*/

void gnuplot_close(gnuplot_ctrl * handle)
{
    int     i ;

    if (pclose(handle->gnucmd) == -1) {
        fprintf(stderr, "problem closing communication to gnuplot\n") ;
        return ;
    }
    if (handle->ntmp) {
        for (i=0 ; i<handle->ntmp ; i++) {
            remove(handle->tmp_filename_tbl[i]) ;
            free(handle->tmp_filename_tbl[i]);
            handle->tmp_filename_tbl[i] = NULL;

        }
    }
    free(handle) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sends a command to an active gnuplot session.
  @param    handle Gnuplot session control handle
  @param    cmd    Command to send, same as a printf statement.

  This sends a string to an active gnuplot session, to be executed.
  There is strictly no way to know if the command has been
  successfully executed or not.
  The command syntax is the same as printf.

  Examples:

  @code
  gnuplot_cmd(g, "plot %d*x", 23.0);
  gnuplot_cmd(g, "plot %.18e * cos(%.18e * x)", 32.0, -3.0);
  @endcode

  Since the communication to the gnuplot process is run through
  a standard Unix pipe, it is only unidirectional. This means that
  it is not possible for this interface to query an error status
  back from gnuplot.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_cmd(gnuplot_ctrl *  handle, char const *  cmd, ...)
{
    va_list ap ;

    va_start(ap, cmd);
    vfprintf(handle->gnucmd, cmd, ap);
    va_end(ap);

    fputs("\n", handle->gnucmd) ;
    fflush(handle->gnucmd) ;
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Change the plotting style of a gnuplot session.
  @param    h Gnuplot session control handle
  @param    plot_style Plotting-style to use (character string)
  @return   void

  The provided plotting style is a character string. It must be one of
  the following:

  - lines
  - points
  - linespoints
  - impulses
  - dots
  - steps
  - errorbars
  - boxes
  - boxeserrorbars
 */
/*--------------------------------------------------------------------------*/

void gnuplot_setstyle(gnuplot_ctrl * h, const char * plot_style)
{
    if (strcmp(plot_style, "lines") &&
        strcmp(plot_style, "points") &&
        strcmp(plot_style, "linespoints") &&
        strcmp(plot_style, "impulses") &&
        strcmp(plot_style, "dots") &&
        strcmp(plot_style, "steps") &&
        strcmp(plot_style, "errorbars") &&
        strcmp(plot_style, "boxes") &&
        strcmp(plot_style, "boxerrorbars")) {
        fprintf(stderr, "warning: unknown requested style: using points\n") ;
        strcpy(h->pstyle, "points") ;
    } else {
        strcpy(h->pstyle, plot_style) ;
    }
    return ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the x label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for X label.
  @return   void

  Sets the x label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_xlabel(gnuplot_ctrl * h, const char * label)
{
    gnuplot_cmd(h, "set xlabel \"%s\"", label) ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Sets the y label of a gnuplot session.
  @param    h Gnuplot session control handle.
  @param    label Character string to use for Y label.
  @return   void

  Sets the y label for a gnuplot session.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_set_ylabel(gnuplot_ctrl * h, const char * label)
{
    gnuplot_cmd(h, "set ylabel \"%s\"", label) ;
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Resets a gnuplot session (next plot will erase previous ones).
  @param    h Gnuplot session control handle.
  @return   void

  Resets a gnuplot session, i.e. the next plot will erase all previous
  ones.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_resetplot(gnuplot_ctrl * h)
{
    int     i ;
    if (h->ntmp) {
        for (i=0 ; i<h->ntmp ; i++) {
            remove(h->tmp_filename_tbl[i]) ;
            free(h->tmp_filename_tbl[i]);
            h->tmp_filename_tbl[i] = NULL;

        }
    }
    h->ntmp = 0 ;
    h->nplots = 0 ;
    return ;
}

/*-------------------------------------------------------------------------*/
/**
  @brief    Plots a 2d graph from a list of doubles.
  @param    handle  Gnuplot session control handle.
  @param    d       Array of doubles.
  @param    n       Number of values in the passed array.
  @param    title   Title of the plot.
  @return   void

  Plots out a 2d graph from a list of doubles. The x-coordinate is the
  index of the double in the list, the y coordinate is the double in
  the list.

  Example:

  @code
    gnuplot_ctrl    *h ;
    double          d[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        d[i] = (double)(i*i) ;
    }
    gnuplot_plot_x(h, d, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_x(
    gnuplot_ctrl    *   handle,
    double          *   d,
    int                 n,
    const char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || d==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n ; i++) {
      fprintf(tmpfd, "%.18e\n", d[i]);
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title);
    return ;
}



/*-------------------------------------------------------------------------*/
/**
  @brief    Plot a 2d graph from a list of points.
  @param    handle      Gnuplot session control handle.
  @param    x           Pointer to a list of x coordinates.
  @param    y           Pointer to a list of y coordinates.
  @param    n           Number of doubles in x (assumed the same as in y).
  @param    title       Title of the plot.
  @return   void

  Plots out a 2d graph from a list of points. Provide points through a list
  of x and a list of y coordinates. Both provided arrays are assumed to
  contain the same number of values.

  @code
    gnuplot_ctrl    *h ;
    double          x[50] ;
    double          y[50] ;
    int             i ;

    h = gnuplot_init() ;
    for (i=0 ; i<50 ; i++) {
        x[i] = (double)(i)/10.0 ;
        y[i] = x[i] * x[i] ;
    }
    gnuplot_plot_xy(h, x, y, 50, "parabola") ;
    sleep(2) ;
    gnuplot_close(h) ;
  @endcode
 */
/*--------------------------------------------------------------------------*/

void gnuplot_plot_xy(
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    const char            *   title
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile(handle,tmpfname,title);
    return ;
}


// Added
void gnuplot_plot_atmpfile_two (gnuplot_ctrl * handle, char const* tmp_filename, char const* title,
				char const* title2)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;
    gnuplot_cmd(handle, "%s \"%s\" index 0 title \"%s\" with %s, \"%s\" index 1 title \"%s\" with %s",
		cmd, tmp_filename, title, handle->pstyle, tmp_filename, title2, handle->pstyle) ;
    handle->nplots++ ;
    return ;
}


void gnuplot_plot_xy_two (
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    int                 n,
    const char            *   title,
    double          *   x2,
    double          *   y2,
    int                 n2,
    const char            *   title2
)
{
    int     i ;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (n<1)|| x2==NULL || y2==NULL || (n2<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x[i], y[i]) ;
    }
    fprintf (tmpfd, "\n\n");
    /* Write data to this file  */
    for (i=0 ; i<n2; i++) {
        fprintf(tmpfd, "%.18e %.18e\n", x2[i], y2[i]) ;
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile_two (handle,tmpfname,title,title2);
}


// Added
void gnuplot_plot_atmpfile_xyz (gnuplot_ctrl * handle, char const* tmp_filename)
{
    char const *    cmd    = "splot";
    gnuplot_cmd(handle, "%s \"%s\" using 1:2:3 notit with %s",
		cmd, tmp_filename, handle->pstyle);
    handle->nplots++ ;
    return ;
}


void gnuplot_plot_xyz (
    gnuplot_ctrl    *   handle,
    double          *   x,
    double          *   y,
    double          *   z,
    int                 nx,
    int                 ny,
    int		            connected
)
{
    int     i, j;
    FILE*   tmpfd ;
    char const * tmpfname;

    if (handle==NULL || x==NULL || y==NULL || (nx<1)|| z==NULL || (ny<1)) return ;

    /* Open temporary file for output   */
    tmpfname = gnuplot_tmpfile(handle);
    tmpfd = fopen(tmpfname, "w");

    if (tmpfd == NULL) {
        fprintf(stderr,"cannot create temporary file: exiting plot") ;
        return ;
    }

    /* Write data to this file  */
    for (i=0 ; i<nx; i++) {
      for (j=0; j < ny; j++)
        fprintf(tmpfd, "%.18e %.18e %.18e\n", x[i], y[j], z [j + ny * i]) ;
      fprintf (tmpfd, "\n");
      if (!connected)
		fprintf (tmpfd, "\n");
    }
    fclose(tmpfd) ;

    gnuplot_plot_atmpfile_xyz (handle,tmpfname);
}


/*-------------------------------------------------------------------------*/
/**
  @brief    Open a new session, plot a signal, close the session.
  @param    title   Plot title
  @param    style   Plot style
  @param    label_x Label for X
  @param    label_y Label for Y
  @param    x       Array of X coordinates
  @param    y       Array of Y coordinates (can be NULL)
  @param    n       Number of values in x and y.
  @return

  This function opens a new gnuplot session, plots the provided
  signal as an X or XY signal depending on a provided y, waits for
  a carriage return on stdin and closes the session.

  It is Ok to provide an empty title, empty style, or empty labels for
  X and Y. Defaults are provided in this case.
 */
/*--------------------------------------------------------------------------*/

void gnuplot_one_function (
  string title,
  string style,
  string label_x,
  string label_y,
  double  *   x,
  double  *   y,
  int         n
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || n<1) return ;

  if ((handle = gnuplot_init()) == NULL) return ;
  if (style!="") {
      gnuplot_setstyle(handle, style.c_str());
  } else {
      gnuplot_setstyle(handle, (char *)"lines");
  }
      gnuplot_set_xlabel(handle, label_x.c_str());
      gnuplot_set_ylabel(handle, label_y.c_str());
  if (y==NULL) {
      gnuplot_plot_x(handle, x, n, title.c_str());
  } else {
      gnuplot_plot_xy(handle, x, y, n, title.c_str());
  }
  
  printf("press 2 x ENTER to continue\n");
  while (getchar()!='\n') {}
  while (getchar()!='\n') {}
  gnuplot_close(handle);
  return ;
}


void gnuplot_one_function_3d (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  string   label_z,
  double  *   x,
  double  *   y,
  void  *    z,
  int         nx,
  int         ny,
  int         connected,
  double      view1,
  double      view2
)
{
  gnuplot_ctrl    *   handle ;

  if ((handle = gnuplot_init()) == NULL) return ;
  if (style!="") {
      gnuplot_setstyle(handle, style.c_str());
  } else {
      gnuplot_setstyle(handle, (char *)"lines");
  }
      gnuplot_set_xlabel(handle, label_x.c_str());
       gnuplot_set_ylabel(handle, label_y.c_str());

  gnuplot_cmd(handle, "set title \"%s\"", title.c_str()) ;
  gnuplot_cmd(handle, "set ticslevel 0") ;
  gnuplot_cmd(handle, "set view %f, %f", view1, view2) ;
  gnuplot_cmd(handle, "set zlabel \"%s\"", label_z.c_str()) ;


  gnuplot_plot_xyz (handle, x, y, (double *) z, nx, ny, connected);

  printf("press 2 x ENTER to continue\n");
  while (getchar()!='\n') {}
  while (getchar()!='\n') {}
  gnuplot_close(handle);
  return ;
}


// Added
void gnuplot_two_functions (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  double  *   x,
  double  *   y,
  int         n,
  string   dataname,
  double  *   x2,
  double  *   y2,
  int         n2,
  string   dataname2 
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || y==NULL || x2==NULL || y2==NULL || n<1 || n2<1)
	return ;

  if ((handle = gnuplot_init()) == NULL)
	return ;

  if (style != "")
	gnuplot_setstyle(handle, style.c_str());
  else
	gnuplot_setstyle(handle, (char *)"lines");

	gnuplot_set_xlabel(handle, label_x.c_str());
    
    gnuplot_set_ylabel(handle, label_y.c_str());
    
  gnuplot_cmd(handle, "set title \"%s\"", title.c_str()) ;

  gnuplot_plot_xy_two (handle, x, y, n, dataname.c_str(), x2, y2, n2, dataname2.c_str());

  printf("press 2 x ENTER to continue\n");
  while (getchar()!='\n') {}
  while (getchar()!='\n') {}

  gnuplot_close(handle);
  return ;
}


// Added:
void gnuplot_one_function_jpg (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  double  *   x,
  double  *   y,
  int         n,
  string   filename
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || n<1) return ;

  if ((handle = gnuplot_init()) == NULL) return ;
  gnuplot_cmd(handle, "set term jpeg") ;
  gnuplot_cmd(handle, "set out \"%s\"", filename.c_str()) ;
 if (style != "") {
      gnuplot_setstyle(handle, style.c_str());
  } else {
      gnuplot_setstyle(handle, (char *)"lines");
  }
      gnuplot_set_xlabel(handle, label_x.c_str());
      gnuplot_set_ylabel(handle, label_y.c_str());
  if (y==NULL) {
      gnuplot_plot_x(handle, x, n, title.c_str());
  } else {
      gnuplot_plot_xy(handle, x, y, n, title.c_str());
  }
  gnuplot_close(handle);
  return ;
}


void gnuplot_one_function_square_jpg (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  double  *   x,
  double  *   y,
  int         n,
  string   filename
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || n<1) return ;

  if ((handle = gnuplot_init()) == NULL) return ;
  gnuplot_cmd(handle, "set term jpeg") ;
  gnuplot_cmd(handle, "set out \"%s\"", filename.c_str()) ;
  gnuplot_cmd(handle, "set size square") ;
 if (style!= "") {
      gnuplot_setstyle(handle, style.c_str());
  } else {
      gnuplot_setstyle(handle, (char *)"lines");
  }
      gnuplot_set_xlabel(handle, label_x.c_str());
      gnuplot_set_ylabel(handle, label_y.c_str());
  if (y==NULL) {
      gnuplot_plot_x(handle, x, n, title.c_str());
  } else {
      gnuplot_plot_xy(handle, x, y, n, title.c_str());
  }
  gnuplot_close(handle);
  return ;
}


void gnuplot_two_functions_jpg (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  double  *   x,
  double  *   y,
  int         n,
  string   dataname,
  double  *   x2,
  double  *   y2,
  int         n2,
  string   dataname2,
  string   filename 
)
{
  gnuplot_ctrl    *   handle ;

  if (x==NULL || y==NULL || x2==NULL || y2==NULL || n<1 || n2<1)
	return ;

  if ((handle = gnuplot_init()) == NULL)
	return ;

  gnuplot_cmd(handle, "set term jpeg") ;
  gnuplot_cmd(handle, "set out \"%s\"", filename.c_str()) ;
  if (style!="")
	gnuplot_setstyle(handle, style.c_str());
  else
	gnuplot_setstyle(handle, (char *)"lines");

	gnuplot_set_xlabel(handle, label_x.c_str());
     gnuplot_set_ylabel(handle, label_y.c_str());
    
  gnuplot_cmd(handle, "set title \"%s\"", title.c_str()) ;

  gnuplot_plot_xy_two (handle, x, y, n, dataname.c_str(), x2, y2, n2, dataname2.c_str());
  gnuplot_close(handle);
  return ;
}


void gnuplot_one_function_3d_jpg (
  string   title,
  string   style,
  string   label_x,
  string   label_y,
  string   label_z,
  double  *   x,
  double  *   y,
  void  *   z,
  int         nx,
  int         ny,
  int	      connected,
  double      view1,
  double      view2,
  string   filename 
)
{
  gnuplot_ctrl    *   handle ;

  if ((handle = gnuplot_init()) == NULL) return ;

  gnuplot_cmd(handle, "set term jpeg") ;
  gnuplot_cmd(handle, "set out \"%s\"", filename.c_str()) ;
  if (style!="") {
      gnuplot_setstyle(handle, style.c_str());
  } else {
      gnuplot_setstyle(handle, (char *)"lines");
  }
      gnuplot_set_xlabel(handle, label_x.c_str());
      gnuplot_set_ylabel(handle, label_y.c_str());

  gnuplot_cmd(handle, "set title \"%s\"", title.c_str()) ;
  gnuplot_cmd(handle, "set ticslevel 0") ;
  gnuplot_cmd(handle, "set view %f, %f", view1, view2) ;
  gnuplot_cmd(handle, "set zlabel \"%s\"", label_z.c_str()) ;


  gnuplot_plot_xyz (handle, x, y, (double *) z, nx, ny, connected);

  gnuplot_close(handle);
  return ;
}



void gnuplot_plot_slope(
    gnuplot_ctrl    *   handle,
    double              a,
    double              b,
    char            *   title
)
{
    char const *    cmd    = (handle->nplots > 0) ? (char *)"replot" : (char *)"plot";
    title                  = (title == NULL)      ? (char *)"(none)" : title;

    gnuplot_cmd(handle, "%s %.18e * x + %.18e title \"%s\" with %s",
                  cmd, a, b, title, handle->pstyle) ;

    handle->nplots++ ;
    return ;
}


void gnuplot_plot_equation(
    gnuplot_ctrl    *   h,
    char            *   equation,
    char            *   title
)
{
    char const *    cmd    = (h->nplots > 0) ? (char *)"replot" : (char *)"plot";
    title                  = (title == NULL)      ? (char *)"(none)" : title;

    gnuplot_cmd(h, (char *)"%s %s title \"%s\" with %s",
                  cmd, equation, title, h->pstyle) ;
    h->nplots++ ;
    return ;
}


int gnuplot_write_x_csv(
    char const * fileName,
    double const * d,
    int n,
    char const * title)
{
    int     i;
    FILE*   fileHandle;

    if (fileName==NULL || d==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e\n", i, d[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_xy_csv(
    char const *        fileName,
    double const    *   x,
    double const    *   y,
    int                 n,
    char const      *   title)
{
    int     i ;
    FILE*   fileHandle;

    if (fileName==NULL || x==NULL || y==NULL || (n<1))
    {
        return -1;
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%.18e, %.18e\n", x[i], y[i]) ;
    }

    fclose(fileHandle) ;

    return 0;
}

int gnuplot_write_multi_csv(
    char const *        fileName,
    double const    **  xListPtr,
    int                 n,
    int                 numColumns,
    char const      *   title)
{
    int     i;
    int     j;
    FILE*   fileHandle;

    if (fileName==NULL || xListPtr==NULL || (n<1) || numColumns <1)
    {
        return -1;
    }

    for (j=0;j<numColumns;j++)
    {
        if (xListPtr[j] == NULL)
        {
            return -1;
        }
    }

    fileHandle = fopen(fileName, "w");

    if (fileHandle == NULL)
    {
        return -1;
    }

    // Write Comment.
    if (title != NULL)
    {
        fprintf(fileHandle, "# %s\n", title) ;
    }

    /* Write data to this file  */
    for (i=0 ; i<n; i++)
    {
        fprintf(fileHandle, "%d, %.18e", i, xListPtr[0][i]) ;
        for (j=1;j<numColumns;j++)
        {
            fprintf(fileHandle, ", %.18e", xListPtr[j][i]) ;
        }
        fprintf(fileHandle, "\n");
    }

    fclose(fileHandle) ;

    return 0;
}

char const * gnuplot_tmpfile(gnuplot_ctrl * handle)
{
    static char const * tmp_filename_template = "gnuplot_tmpdatafile_XXXXXX";
    char *              tmp_filename = NULL;
    int                 tmp_filelen = strlen(tmp_filename_template);

#ifndef WIN32
    int                 unx_fd;
#endif // #ifndef WIN32

    assert(handle->tmp_filename_tbl[handle->ntmp] == NULL);

    /* Open one more temporary file? */
    if (handle->ntmp == GP_MAX_TMP_FILES - 1) {
        fprintf(stderr,
                "maximum # of temporary files reached (%d): cannot open more",
                GP_MAX_TMP_FILES) ;
        return NULL;
    }

    tmp_filename = (char*) malloc(tmp_filelen+1);
    if (tmp_filename == NULL)
    {
        return NULL;
    }
    strcpy(tmp_filename, tmp_filename_template);

#ifdef WIN32
    if (_mktemp(tmp_filename) == NULL)
    {
        return NULL;
    }
#else // #ifdef WIN32
    unx_fd = mkstemp(tmp_filename);
    if (unx_fd == -1)
    {
        return NULL;
    }
    close(unx_fd);

#endif // #ifdef WIN32

    handle->tmp_filename_tbl[handle->ntmp] = tmp_filename;
    handle->ntmp ++;
    return tmp_filename;
}

void gnuplot_plot_atmpfile(gnuplot_ctrl * handle, char const* tmp_filename, char const* title)
{
    char const *    cmd    = (handle->nplots > 0) ? "replot" : "plot";
    title                  = (title == NULL)      ? "(none)" : title;
    gnuplot_cmd(handle, "%s \"%s\" title \"%s\" with %s", cmd, tmp_filename,
                  title, handle->pstyle) ;
    handle->nplots++ ;
    return ;
}


/* vim: set ts=4 et sw=4 tw=75 */
