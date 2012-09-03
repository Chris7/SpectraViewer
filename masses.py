"""
Author: Chris Mitchell (chris.mit7@gmail.com)
Copyright (C) 2012 Chris Mitchell

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
#masses from: http://www.weddslist.com/ms/tables.html#tm4
#protein weights with some mascot custom mods too, stored as weight,charge
protein_weights =  {'G': (57.021464, 0),
                    'A': (71.037114, 0),
                    'S': (87.032028, 0),
                    'P': (97.052764, 0),
                    'V': (99.068414, 0),
                    'T': (101.047678, 0),
                    'C': (103.009184, 0),
                    'I': (113.084064, 0),
                    'L': (113.084064, 0),
                    'N': (114.042927, 0),
                    'D': (115.026943, 0),
                    'Q': (128.058578, 0),
                    'K': (128.094963, 1),
                    'E': (129.042593, 0),
                    'M': (131.040485, 0),
                    'H': (137.058912, 1),
                    'F': (147.068414, 0),
                    'R': (156.101111, 1),
                    'Y': (163.063329, 0),
                    'W': (186.079313, 0)
                    }

mod_weights = {'h': (1.007825, 0),
               'h2o': (18.010565, 0),
               'h2o': (18.0106, 0),
               'nh3': (17.026549, 0),
               'ch2': (14.015650, 0),
               'methylation': (14.015650, 0),
               'o': (15.994915, 0),
               'oxidation': (15.994915, 0),
               'acetylation': (42.010565, 0),
               'acetylation': (42.0106, 0),
               'carbamidation': (57.021464, 0),
               'carboxylation': (58.005479, 0),
               'phosphorylation': (79.966330, 0),
               'amidation': (0.984016, 0),
               'formylation': (27.994915, 0),
               'cho': (29.002739665, 0),
               'nh2': (16.01872407, 0),
               'co': (27.99491463, 0),
               'oh': (17.00274, 0)
               }

lossMasses = {'K': ((-1*mod_weights['nh3'][0],('a','b','y'),'-nh3'),(mod_weights['h2o'][0],('b',),'+h2o')),
              'R': ((-1*mod_weights['nh3'][0],('a','b','y'),'-nh3'),(mod_weights['h2o'][0],('b',),'+h2o')),
              'S': ((-1*mod_weights['h2o'][0],('b'),'-h2o'),),
              'T': ((-1*mod_weights['h2o'][0],('b'),'-h2o'),),
              'D': ((-1*mod_weights['h2o'][0],('b','y'),'-h2o'),),
              'E': ((-1*mod_weights['h2o'][0],('b','y'),'-h2o'),)}