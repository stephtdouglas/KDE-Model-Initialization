import glob, os
import astropy.io.ascii as at
import numpy as np
from astropy.table import vstack


if __name__=="__main__":
    # First, combine Angie's individual model files into one
    models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/UpSco_KDE_init/*.txt"))
    models = np.sort(models)

    model_tables = []
    for mfile in models:
        temp_tab = at.read(mfile)
        mass = float(mfile.split("/")[-1][:3])
        temp_tab["Mass"] = np.ones(len(temp_tab))*mass
        model_tables.append(temp_tab)
    mod = vstack(model_tables)
    mod.rename_column("Period","prot")

    mod.write("all_KDE.csv",delimiter=",",overwrite=True)

    # Now, need to account for the stars where Sean changed the initial periods
    # because they were rotating faster than breakup
    model = "WideHat8Myr_Mattea2015"
    model_age = 8
    mod_file = os.path.expanduser(f"~/Dropbox/Models/{model}/{model}_{model_age:05d}Myr.txt")
    mod2 = at.read(mod_file,names=["mass","prot"])
    tab = mod.copy()
    tab["Sean_prot"] = mod2["prot"]
    tab["Sean_mass"] = mod2["mass"]
    tab["pdiff"] = tab["prot"] - tab["Sean_prot"]
    tab["mdiff"] = tab["Mass"] - tab["Sean_mass"]
    tab["pdiff_frac"] = np.abs(tab["prot"] - tab["Sean_prot"])/tab["prot"]
    tab["mdiff_frac"] = np.abs(tab["Mass"] - tab["Sean_mass"])/tab["Mass"]

    # Where the periods differ significantly, set probability to 0 (so they
    # won't be included in anything), and redistribute that probability to
    # other stars in the same mass bin
    mass_bins = np.arange(0.05,1.4,0.1)
    for i in range(len(mass_bins)-1):
        sub = tab[500*i:500*(i+1)]
        min_per = min(sub["Sean_prot"])
        # This is marking the wrong things as bad. AAAAAAAAAAAAHHHH
        # bad = sub["pdiff"]<0.005
        bad = sub["pdiff_frac"]>0.05
        nbad = len(np.where(bad)[0])
        pbad = np.sum(sub["Prob"][bad])
        # print(mass_bins[i],nbad,pbad)
        ngood = 500-nbad

        # print(np.sum(sub["Prob"]))

        sub["Prob"][~bad] += (pbad/ngood)
        sub["Prob"][bad] = 0

        # print(np.sum(sub["Prob"]))

    tab.write("all_KDE_corrected.csv",delimiter=",",overwrite=True)
