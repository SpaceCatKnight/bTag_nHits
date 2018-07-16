import FinalClusterMatcher as FCM
import sys


if __name__ == "__main__":

	if len(sys.argv) > 1: i = sys.argv[1]

	BG_noPU_String = 'root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/store/user/thaarres/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/btagHits_wPVs/180502_130824/0000/flatTuple_{}.root' #1-47

	FCM.dR_Dist('BG_{}'.format(i),[BG_noPU_String.format(i)], 350, BG=True, EarlyBreak=0)

