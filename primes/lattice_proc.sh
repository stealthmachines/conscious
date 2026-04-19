#!/bin/sh
# lattice_proc.sh -- phi4096 lattice-native process scheduler
# Seed: 0xfc3ed6b03ab19179   Slots: 4096

echo '[lattice-proc] Installing scheduler tools...'
apk add --quiet util-linux schedutils 2>/dev/null || apk add --quiet util-linux 2>/dev/null

echo '[lattice-proc] Configuring cgroup v2 lattice_procs...'
mount | grep -q cgroup2 || mount -t cgroup2 none /sys/fs/cgroup 2>/dev/null
mkdir -p /sys/fs/cgroup/lattice_procs 2>/dev/null
echo '+cpu +memory +io' > /sys/fs/cgroup/cgroup.subtree_control 2>/dev/null
echo '1' > /sys/fs/cgroup/lattice_procs/cpu.weight 2>/dev/null
echo '25000 100000' > /sys/fs/cgroup/lattice_procs/cpu.max 2>/dev/null
echo '64M' > /sys/fs/cgroup/lattice_procs/memory.max 2>/dev/null

echo '[lattice-proc] Writing /usr/local/bin/sched_run...'
cat > /usr/local/bin/sched_run << 'SCHED_RUN_EOF'
#!/bin/sh
# lattice-native process launcher -- spawns $@ under lattice scheduler
exec chrt -o 0 "$@"
SCHED_RUN_EOF
chmod +x /usr/local/bin/sched_run

echo '[lattice-proc] Writing /etc/profile.d/lattice_sched.sh...'
cat > /etc/profile.d/lattice_sched.sh << 'PROF_EOF'
export SCHED_POLICY=SCHED_OTHER
export SCHED_NICE=-20
export SCHED_AFFINITY=0x01
export IONICE_CLASS=none
export IONICE_LEVEL=0
export CGROUP_WEIGHT=1
export MEM_LIMIT=64M
alias run='sched_run'
echo "[lattice-sched] policy=$SCHED_POLICY  affinity=$SCHED_AFFINITY  cgroup-weight=$CGROUP_WEIGHT"
PROF_EOF

echo '[lattice-proc] Applying scheduler to current shell (PID $$)...'
chrt -o -p 0 $$ 2>/dev/null && echo '  policy: SCHED_OTHER' || echo '  [warn] chrt unavailable'
taskset -p 0x01 $$ 2>/dev/null && echo '  affinity: 0x01' || echo '  [warn] taskset unavailable'
renice -n -20 -p $$ 2>/dev/null && echo '  nice: -20' || true
echo '-500' > /proc/$$/oom_score_adj 2>/dev/null && echo '  oom_score_adj: -500' || true
echo $$ > /sys/fs/cgroup/lattice_procs/cgroup.procs 2>/dev/null && echo '  cgroup: lattice_procs (weight=1)' || true
ulimit -s 64 2>/dev/null && echo '  stack: 64 kB' || true

echo '[lattice-proc] Lattice-native process environment active.'
echo "  Run any command under the lattice scheduler: sched_run <cmd>"
echo "  Or use the alias:  run <cmd>"
