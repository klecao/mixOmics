# Copyright (C) 2009 
# S�bastien D�jean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonz�lez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh L� Cao, French National Institute for Agricultural Research and 
# Queensland Facility for Advanced Bioinformatics, University of Queensland, Australia
# Pierre Monget, Ecole d'Ingenieur du CESI, Angouleme, France
#
# This function was borrowed from the mclust package
#  
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


map <-
function (Y, ...) 
{
    nrowY <- nrow(Y)
    cl <- numeric(nrowY)
    I <- 1:nrowY
    J <- 1:ncol(Y)
    for (i in I) {
        cl[i] <- (J[Y[i, ] == max(Y[i, ])])[1]
    }
    cl
}