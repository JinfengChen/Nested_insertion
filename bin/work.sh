echo "HEG4 double pong reads"
grep "14729:99648" all.blatout | less -S
echo "identify double pong/ping"
python ReNameSRA_RelocaTEi_Nest_insertion.py --input Run_folder

