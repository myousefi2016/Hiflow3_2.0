#!/usr/bin/perl -w

# Copyright (C) 2011-2017 Vincent Heuveline
#
# HiFlow3 is free software: you can redistribute it and/or modify it under the
# terms of the European Union Public Licence (EUPL) v1.2 as published by the
#/ European Union or (at your option) any later version.
#
# HiFlow3 is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the European Union Public Licence (EUPL) v1.2 for more
# details.
#
# You should have received a copy of the European Union Public Licence (EUPL) v1.2
# along with HiFlow3.  If not, see <https://joinup.ec.europa.eu/page/eupl-text-11-12>.

# \brief Converts the LOG_INFO file to an HTML overview
# \author Thomas Gengenbach

######################################################################
############################## MAIN ##################################
######################################################################

use strict;

unless (@ARGV eq 3) {
    print "Usage: info2html.pl <log-file> <path-to-hiflow3> <application-name>\n\n";
    die("Wrong number of arguments!\n");
}

my $APPLICATION = "$ARGV[$#ARGV]";
print "HTML data sheet is generated for $APPLICATION.\n";

# Set path to hiflow info.log file
my $LOG_INFO_FILE = "$ARGV[0]";

my @path_to_info_log = "";
# if name contains path
if ($LOG_INFO_FILE =~ /\//){
    print "The string contains non-alphanumeric characters";
    @path_to_info_log = split(/\//, $LOG_INFO_FILE);
    pop(@path_to_info_log);
} else {
    @path_to_info_log = ".";
}

my $path_to_info_log = "";
foreach (@path_to_info_log) {
    $path_to_info_log .=  $_ . "/";
}
chop($path_to_info_log);

unless (-e "$LOG_INFO_FILE") {
    die("Log file $LOG_INFO_FILE doesn:t exist!");
}

# Set path to HiFlow3 repository
my $path_to_hiflow = "$ARGV[1]";

my $APP = "";
my @app_files;
application();

my $MESH = "";
my @mesh_files;
mesh();

my $DOF = "";
my @dof_files;
dof();

my $FEM = "";
my @fem_files;
fem();

my $LA = "";
my @la_files;
linear_algebra();

my $LS = "";
my @ls_files;
linear_solver();

my $NL = "";
my @nl_files;
nonlinear_solver();

my $VS = "";
my $VISU = "";
my @vs_files;
my @visu_files;
vector_space();

# Parsing
my $INFO = "";
my $date = "";
my $i = 0;

my $plot_it = 0;
parse_info_log();

my $residual = "";
if ($plot_it) {
    plot("Residual");
}

if ($VISU eq "") {
    $VISU .= "<tr><td>Visualization</td>";
    $VISU .= "<td>Either none or not with the Visualization class</td></tr>";
}
write_html();

print "Written file $path_to_info_log/hiflow3_$APPLICATION.html.\n";
print "Exit succesfully!\n";

exit 0;

######################################################################
########################### SUBROUTINES ##############################
######################################################################

######################################################################

# write HTML
sub write_html {

######################################################################
## Write HTML file
    open(HTML, ">$path_to_info_log/hiflow3_$APPLICATION.html");
    print HTML <<ENDOFHTML;
    <html>
	<head>
	<title>[HiFlow3 - Simulation Report]</title>
	</head>

	<body>
	<center><image src="$path_to_hiflow/utils/hiflow3_logo.png" style="width:30%"></center>

	<!-- ------------------------------------------------------------- -->
	<!-- CONTENTS                                                      -->
	<!-- ------------------------------------------------------------- -->

       	<h1><a href="http://www.hiflow3.org">HiFlow3</a> Application $APPLICATION</h1>

	<table border=1 cellspacing=0 cellpadding=2>
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>Contents of Report</b></td>
	</tr>
	<tr>
	<td><p align=right>I.</p></td>
	<td><a href="#Parameters">Parameters</a></td>
	</tr>
	<tr>
	<td><p align=right>II.</p></td>
	<td><a href="#Modules">Modules</a></td>
	</tr>
	<tr>
	<td><p align=right>III.</p></td>
	<td><a href="#Visualization">Visualization</a></td>
	</tr>
	<tr>
	<td><p align=right>IV.</p></td>
	<td><a href="#Summary">Summary</a></td>
	</tr>
	</table>

	<!-- ------------------------------------------------------------- -->
	<br><div style='border:none;border-bottom:solid windowtext .75pt;width:100%'></div>

	<!-- ------------------------------------------------------------- -->
	<!-- PARAMETERS                                                    -->
	<!-- ------------------------------------------------------------- -->
	<h2><a name=Parameters></a>I. Parameters</h2>

	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>Meta Data</b></td>
	</tr>
	<tr>
	<td><p>Version</p></td>
	<td><p>0.1</p></td>
	</tr>
	<tr>
	<td><p>Simulation date</p></td>
	<td><p>$date</p></td>
	</tr>
	<tr>
	<td><p>Created from log file</p></td>
	<td><p>$LOG_INFO_FILE</p></td>
	</tr>
	</table>
	<br><br>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>Simulation Parameters</b></td>
	</tr>
	$APP
	</table>

	<br>
	<a href="#top">top</a><br>
	<!-- ------------------------------------------------------------- -->
	<br><div style='border:none;border-bottom:solid windowtext .75pt;width:100%'></div>

	<!-- ------------------------------------------------------------- -->
	<!-- MODULES                                                       -->
	<!-- ------------------------------------------------------------- -->

	<h2><a name=Modules></a>II. Modules</h2>
ENDOFHTML
    my $module_count = 0;
    if ($MESH ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;

	<!-- ------------------------------------------------------------- -->
	<!-- MESH                                                          -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Mesh></a>II.$module_count Mesh</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>Mesh Parameters</b></td>
	</tr>
	$MESH
	</table>
	<br>
ENDOFHTML
}
    if ($DOF ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;

	<!-- ------------------------------------------------------------- -->
	<!-- DOF                                                           -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Degree of Freedom></a>II.$module_count Degree of Freedom</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>DoF Parameters</b></td>
	</tr>
	$DOF
	</table>
	<br>
ENDOFHTML
    }
    if ($FEM ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;

	<!-- ------------------------------------------------------------- -->
	<!-- FEM                                                           -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Finite Elements></a>II.$module_count Finite Elements</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>FEM Parameters</b></td>
	</tr>
	$FEM
	</table>
	<br>
ENDOFHTML
    }
    if ($LA ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;

	<!-- ------------------------------------------------------------- -->
	<!-- Linear Algebra                                                -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Linear Algebra></a>II.$module_count Linear Algebra</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>LA Parameters</b></td>
	</tr>
	$LA
	</table>
	<br>
ENDOFHTML
    }
    if ($LS ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;

	<!-- ------------------------------------------------------------- -->
	<!-- Linear Solver                                                 -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Linear Solver></a>II.$module_count Linear Solver</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>LS Parameters</b></td>
	</tr>
	$LS
	</table>
ENDOFHTML
    }
    if ($NL ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;
	<!-- ------------------------------------------------------------- -->
	<!-- Nonlinear Solver                                              -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Nonlinear Solver></a>II.$module_count Nonlinear Solver</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>NLS Parameters</b></td>
	</tr>
	$NL
	</table>
	<br>
ENDOFHTML
    }
    if ($VS ne "") {
	++$module_count;
        print HTML <<ENDOFHTML;
	<!-- ------------------------------------------------------------- -->
	<!-- Nonlinear Solver                                              -->
	<!-- ------------------------------------------------------------- -->

	<h3><a name=Nonlinear Solver></a>II.$module_count Vector Space</h3>
	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>VS Parameters</b></td>
	</tr>
	$VS
	</table>
	<br>
ENDOFHTML
    }
        print HTML <<ENDOFHTML;

	<a href="#top">top</a><br>
	<!-- ------------------------------------------------------------- -->
	<br><div style='border:none;border-bottom:solid windowtext .75pt;width:100%'></div>

	<!-- ------------------------------------------------------------- -->
	<!-- VISUALIZATION                                                 -->
	<!-- ------------------------------------------------------------- -->

	<h2><a name=Visualization></a>III. Visualization</h2>

	<table border=1 cellspacing=0 cellpadding=2 width=100%>
	<colgroup width=30% />
	<tr>
	<td colspan=2 style='background:#D9D9D9'><b>Visualization Parameters</b></td>
	</tr>
	$VISU
	</table>
	<br>

	<a href="#top">top</a><br>
	<!-- ------------------------------------------------------------- -->
	<br><div style='border:none;border-bottom:solid windowtext .75pt;width:100%'></div>

	<!-- ------------------------------------------------------------- -->
	<!-- SUMMARY                                                       -->
	<!-- ------------------------------------------------------------- -->

	<h2><a name=Summary></a>IV. Performance Summary</h2>
	<br>
ENDOFHTML
    if ($residual ne "") {
	if (-e "$path_to_info_log/Residual.png") {
	    print HTML <<ENDOFHTML;
	    <h3>Residual plot</h3>
	    <img src="$path_to_info_log/Residual.png" alt='Residual Plot'/>
ENDOFHTML
	} else {
	    print HTML <<ENDOFHTML;
	    <h3>Residual plot</h3>
	    <PRE>$residual</PRE>
ENDOFHTML
	}
    }
    print HTML <<ENDOFHTML;
	<a href="#top">top</a><br>
	<!-- ------------------------------------------------------------- -->
	<br><div style='border:none;border-bottom:solid windowtext .75pt;width:100%'></div>

	<!-- ------------------------------------------------------------- -->
	<!-- APPENDIX                                                      -->
	<!-- ------------------------------------------------------------- -->
	<br>
	<center>This HiFlow3 report was created with the Perl script info2html.pl.
	<br>
	NumHPC/EMCL, July 2010.
	</center>

	</body>

	</html>
ENDOFHTML
    close(HTML);
}

######################################################################

######################################################################

# \brief Parse info log file written by HiFlow3.
#
# Needs the filenames from the CMakeLists.txt files.
sub parse_info_log {

######################################################################
    my $residual_tag = "Residual";
    if (-e "$path_to_info_log/$residual_tag.dat") {
	die("File $path_to_info_log/$residual_tag.dat already exists. Rename or delete it!\n");
    }
    $residual_tag .= ":";

    open(INFO, $LOG_INFO_FILE);

    while (<INFO>) {
	my @current_line = split(/ /, $_);
	my @filename = split(/\./, $current_line[0]);

 	my @label_split = split(/:/, $_);
 	my @label;
	my @data;
	if (@label_split ge 3) {
	    @label = split(/ /, $label_split[1]);
	    @data = split(/ /, $label_split[2]);

	    shift(@label);
	    shift(@data);
	}
	my $filename = substr($filename[0], 1);

	if(/INFO/) {
	    $INFO .= "@label";
	    $INFO .= "@data";
	    $INFO .= "<br>";
	    shift(@current_line);
	    shift(@current_line);
	    shift(@current_line);
	    $date .= "@current_line";
	}
	if (grep { $_ eq $filename} @app_files) {
	    $APP .= "<tr><td>@label</td>";
	    $APP .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @mesh_files) {
	    $MESH .= "<tr><td>@label</td>";
	    $MESH .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @dof_files) {
	    $DOF .= "<tr><td>@label</td>";
	    $DOF .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @fem_files) {
	    $FEM .= "<tr><td>@label</td>";
	    $FEM .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @la_files) {
	    $LA .= "<tr><td>@label</td>";
	    $LA .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @ls_files) {
	    if ($current_line[1] eq $residual_tag) {
		$plot_it = 1;
		$filename = $current_line[1];
		chop($filename);
		dat($filename, $current_line[2]);
	    } else {
		$LS .= "<tr><td>@label</td>";
		$LS .= "<td>@data</td></tr>";
	    }
	}
	if (grep { $_ eq $filename} @nl_files) {
	    $NL .= "<tr><td>@label</td>";
	    $NL .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @vs_files) {
	    $VS .= "<tr><td>@label</td>";
	    $VS .= "<td>@data</td></tr>";
	}
	if (grep { $_ eq $filename} @visu_files) {
	    $VISU .= "<tr><td>@label</td>";
	    $VISU .= "<td>@data</td></tr>";
	}
    }
    close(INFO);
}

######################################################################

# create .dat file for GNUPLOT
sub dat {

######################################################################
    open(GP, ">>$path_to_info_log/$_[0].dat");
    print GP "$i $_[1]";
    close(GP);
    ++$i;
}

######################################################################

# create GNUPLOT graphic
sub plot {

######################################################################
    if (-e "$path_to_info_log/$_[0].txt") {
	die("File $path_to_info_log/$_[0].txt already exists. Rename or delete it!\n");
    }

    # ASCII
    open (GNUPLOT, "|gnuplot");
    print GNUPLOT <<EOPLOT;
    set term dumb
	set output "$path_to_info_log/$_[0].txt"
	set data style line
	set xlabel "Number of Iterations"
	set ylabel "Residual"
	set title "$_[0]"

	plot "$path_to_info_log/$_[0].dat" using 1:2 w lines 1
EOPLOT
    close(GNUPLOT);

    # PNG
    open (GNUPLOT, "|gnuplot");
    print GNUPLOT <<EOPLOT;
    set term png
	set output "$path_to_info_log/$_[0].png"
	set data style line
	set xlabel "Number of Iterations"
	set ylabel "Residual"
	set title "$_[0]"

	plot "$path_to_info_log/$_[0].dat" using 1:2 w lines 1
EOPLOT
    close(GNUPLOT);

    system("rm $path_to_info_log/$_[0].dat");

    $residual = "";
    open(RES, "$path_to_info_log/Residual.txt");
    while(<RES>) {
	$residual .= $_;
    }
    close(RES);
    system("rm $path_to_info_log/Residual.txt");
}

######################################################################

# Gather application files to compare names in info.log
sub application {

######################################################################
    $APP = "";
    my $dir = "$path_to_hiflow/test/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@app_files, $src[0]);
	}
    }
}

######################################################################

# Gather mesh files to compare names in info.log
sub mesh {

######################################################################
    my $dir = "$path_to_hiflow/src/mesh/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@mesh_files, $src[0]);
	}
    }
}

######################################################################

# Gather DoF files to compare names in info.log
sub dof {

######################################################################
    my $dir = "$path_to_hiflow/src/dof/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@dof_files, $src[0]);
	}
    }
}

######################################################################

# Gather fem files to compare names in info.log
sub fem {

######################################################################
    my $dir = "$path_to_hiflow/src/fem/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@fem_files, $src[0]);
	}
    }
}

######################################################################

# Gather linear algebra files to compare names in info.log
sub linear_algebra {

######################################################################
    my $dir = "$path_to_hiflow/src/linear_algebra/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@la_files, $src[0]);
	}
    }
}

######################################################################

# Gather linear solver files to compare names in info.log
sub linear_solver {

######################################################################
    my $dir = "$path_to_hiflow/src/linear_solver/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@ls_files, $src[0]);
	}
    }
}

######################################################################

# Gather nonlinear solver files to compare names in info.log
sub nonlinear_solver {

######################################################################
    my $dir = "$path_to_hiflow/src/nonlinear/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    push(@nl_files, $src[0]);
	}
    }
}

######################################################################

# Gather vector space files to compare names in info.log
sub vector_space {

######################################################################
    my $dir = "$path_to_hiflow/src/space/";

    opendir DIR, $dir or die "cannot open dir $dir: $!";
    my @src_names = readdir DIR;
    closedir DIR;

    foreach (@src_names) {
	if ($_ ne "." and $_ ne "..") {
	    my @src = split(/\./, $_);
	    if ($src[0] eq "visualization") {
		push(@visu_files, $src[0]);
	    } else {
		push(@vs_files, $src[0]);
	    }
	}
    }
}
