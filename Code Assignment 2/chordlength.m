function [chord] = chordlength(span,rootchord,tipchord,radial_coord);
    chord = rootchord - (rootchord - tipchord)*(radial_coord/span);
end