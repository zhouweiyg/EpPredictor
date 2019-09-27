#nohup python -u ep_lower.py --cellline 'IMR90' --cols '20' >IMR90_lower.log 2>&1 &
nohup python -u ep_lower.py --cellline 'GM12878' --cols '20' >GM12878_lower.log 2>&1 &
nohup python -u ep_lower.py --cellline 'K562' --cols '20' >K562_lower.log 2>&1 &

#nohup python -u ep_lower.py --cellline 'Hela-S3' --cols '20' >Hela-S3_lower.log 2>&1 &
#nohup python -u ep_lower.py --cellline 'HUVEC' --cols '20' >HUVEC_lower.log 2>&1 &
#nohup python -u ep_lower.py --cellline 'NHEK' --cols '20' >NHEK_lower.log 2>&1 &