# Colab 初始化代码（每次换运行时后运行一次）。运行时自带的 sshd 监听在 127.0.0.1:2222，
# 这里只加入 claude-colab 公钥（sshd 只允许密钥登录），再用 cloudflared 临时隧道把 2222 暴露出来。
# 运行完把打印出的 HOST 发给 Claude（它会写入 research/colab/HOST）。
PUBKEY = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIEJYLo2MiNRdxMKUS3YuG0YujDzgmdXnSHkN1YtFvIKC claude-colab"
import subprocess, os, re, time, socket
def sh(c): subprocess.run(c, shell=True, check=True)
os.makedirs("/root/.ssh", exist_ok=True); os.chmod("/root/.ssh", 0o700)
ak = "/root/.ssh/authorized_keys"
if PUBKEY not in (open(ak).read() if os.path.exists(ak) else ""):
    with open(ak, "a") as f: f.write(PUBKEY + "\n")
os.chmod(ak, 0o600)
def up(port):
    try: socket.create_connection(("127.0.0.1", port), 2).close(); return True
    except OSError: return False
if not up(2222):
    sh("mkdir -p /run/sshd && (ls /etc/ssh/ssh_host_*_key >/dev/null 2>&1 || ssh-keygen -A) && /usr/sbin/sshd")
if not os.path.exists("/usr/local/bin/cloudflared"):
    sh("wget -q -O /usr/local/bin/cloudflared https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64 && chmod +x /usr/local/bin/cloudflared")
subprocess.run("pkill cloudflared", shell=True); time.sleep(1)
subprocess.Popen("nohup cloudflared tunnel --no-autoupdate --url ssh://127.0.0.1:2222 > /content/cf.log 2>&1 &", shell=True)
m = None
for _ in range(90):
    time.sleep(1)
    if os.path.exists("/content/cf.log"):
        m = re.search(r"https://([a-z0-9-]+\.trycloudflare\.com)", open("/content/cf.log").read())
        if m: break
print("HOST:", m.group(1) if m else "not found, see /content/cf.log")
print("CPUs:", os.cpu_count(), " RAM GB:", round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2**30, 1))
