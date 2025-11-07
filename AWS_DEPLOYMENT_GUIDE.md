# AWS EC2 Deployment Guide with Conda

Complete step-by-step guide for deploying the Enrich_KEs Streamlit app to AWS EC2 using Conda.

---

## Prerequisites

Before starting, make sure you have:

- ✅ AWS EC2 instance running (Amazon Linux 2 or similar)
- ✅ SSH access to your EC2 instance
- ✅ Security key (.pem file)
- ✅ GitHub repository set up (https://github.com/mei-franzi/Enrich_KEs)

---

## Part 1: Initial EC2 Setup

### Step 1.1: Connect to EC2

```bash
# From your local machine
ssh -i ~/path/to/your-key.pem ec2-user@your-ec2-ip-address
```

### Step 1.2: Update System

```bash
# Update all packages
sudo yum update -y

# Install essential tools
sudo yum install git wget -y
```

---

## Part 2: Install Miniconda

### Step 2.1: Download Miniconda

```bash
# Navigate to home directory
cd ~

# Download Miniconda installer for Linux
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Verify download (optional)
ls -lh Miniconda3-latest-Linux-x86_64.sh
```

### Step 2.2: Install Miniconda

```bash
# Run installer (silent mode, installs to ~/miniconda3)
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

# Initialize conda for bash shell
~/miniconda3/bin/conda init bash

# Reload shell configuration
source ~/.bashrc

# Verify conda installation
conda --version
# Should output: conda 24.x.x or similar

# Update conda
conda update -n base -c defaults conda -y
```

### Step 2.3: Configure Conda (Optional but Recommended)

```bash
# Disable auto-activation of base environment
conda config --set auto_activate_base false

# Add conda-forge and bioconda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Verify channels
conda config --show channels
```

---

## Part 3: Clone Repository and Setup Environment

### Step 3.1: Clone Your Repository

```bash
# Navigate to home directory
cd ~

# Clone your GitHub repository
git clone https://github.com/mei-franzi/Enrich_KEs.git

# Navigate into the project
cd Enrich_KEs

# Verify files are present
ls -la
```

### Step 3.2: Create environment.yml File

```bash
# Create the environment configuration file
nano environment.yml
```

**Paste this content:**

```yaml
name: enrich_kes
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - pip
  
  # Future bioinformatics tools can be added here:
  # Examples (uncomment when needed):
  # - fastqc
  # - multiqc
  # - star
  # - samtools
  # - bedtools
  # - deeptools
  
  - pip:
      - streamlit>=1.28.0
      - pandas>=2.0.0
      - numpy>=1.24.0
      - scipy>=1.10.0
      - statsmodels>=0.14.0
      - openpyxl>=3.1.0
      - plotly>=5.14.0
      - gprofiler-official>=1.0.0
      - matplotlib>=3.7.0
      - pytest>=7.4.0
      - pytest-cov>=4.1.0
```

**Save and exit:** Press `Ctrl+X`, then `Y`, then `Enter`

### Step 3.3: Create Conda Environment

```bash
# Create environment from yml file
conda env create -f environment.yml

# This will take 5-10 minutes - be patient!
```

### Step 3.4: Activate Environment

```bash
# Activate the environment
conda activate enrich_kes

# Verify activation (prompt should show: (enrich_kes))
which python
# Should output: /home/ec2-user/miniconda3/envs/enrich_kes/bin/python

# List installed packages
conda list

# Test Python
python --version
# Should output: Python 3.10.x
```

### Step 3.5: Verify Installation

```bash
# Test imports
python -c "import streamlit; import pandas; import plotly; print('✅ All imports successful!')"

# Test Streamlit
streamlit --version
```

---

## Part 4: Configure AWS Security Group

**You need to allow traffic on port 8501 (Streamlit's port)**

### Step 4.1: Via AWS Console

1. Go to **AWS Console** → **EC2** → **Instances**
2. Click on your instance
3. Click on **Security** tab
4. Click on the **Security Group** link
5. Click **Edit inbound rules**
6. Click **Add rule**
7. Configure:
   - **Type**: Custom TCP
   - **Port range**: 8501
   - **Source**: 
     - For testing: "My IP" (your current IP)
     - For production: "Anywhere" (0.0.0.0/0)
   - **Description**: Streamlit app
8. Click **Save rules**

### Step 4.2: Verify Security Group (Optional)

```bash
# From your EC2 terminal, check if port is accessible
sudo yum install nmap -y
nmap -p 8501 localhost
```

---

## Part 5: Test Run Your App

### Step 5.1: Manual Test Run

```bash
# Make sure you're in the project directory
cd ~/Enrich_KEs

# Activate conda environment
conda activate enrich_kes

# Run Streamlit
streamlit run app.py --server.port 8501 --server.address 0.0.0.0

# You should see output like:
#   You can now view your Streamlit app in your browser.
#   Network URL: http://0.0.0.0:8501
#   External URL: http://your-ip:8501
```

### Step 5.2: Access Your App

1. **Get your EC2 public IP:**
   ```bash
   curl http://checkip.amazonaws.com
   ```

2. **Open in browser:**
   ```
   http://your-ec2-public-ip:8501
   ```

3. **Test the app:**
   - Upload a file
   - Run enrichment
   - Check all features work

4. **Stop the test run:**
   - Press `Ctrl+C` in the terminal

---

## 🔄 Part 6: Production Setup (Keep App Running)

### Option A: Using Screen (Simple, Good for Testing)

```bash
# Install screen
sudo yum install screen -y

# Create a new screen session
screen -S streamlit_app

# Activate conda and run app
conda activate enrich_kes
cd ~/Enrich_KEs
streamlit run app.py --server.port 8501 --server.address 0.0.0.0

# Detach from screen: Press Ctrl+A, then D
# Your app keeps running!

# Later, to reattach:
screen -r streamlit_app

# To list all screens:
screen -ls

# To kill the screen:
screen -X -S streamlit_app quit
```

### Option B: Using systemd (Recommended for Production)

#### Step 6.1: Create Systemd Service File

```bash
# Create service file
sudo nano /etc/systemd/system/streamlit.service
```

**Paste this content:**

```ini
[Unit]
Description=Streamlit Enrich_KEs App
After=network.target

[Service]
Type=simple
User=ec2-user
WorkingDirectory=/home/ec2-user/Enrich_KEs
Environment="PATH=/home/ec2-user/miniconda3/envs/enrich_kes/bin:/home/ec2-user/miniconda3/bin:/usr/local/bin:/usr/bin:/bin"
ExecStart=/home/ec2-user/miniconda3/envs/enrich_kes/bin/streamlit run app.py --server.port 8501 --server.address 0.0.0.0
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

**Save and exit:** Press `Ctrl+X`, then `Y`, then `Enter`

#### Step 6.2: Enable and Start Service

```bash
# Reload systemd to recognize new service
sudo systemctl daemon-reload

# Enable service to start on boot
sudo systemctl enable streamlit

# Start the service
sudo systemctl start streamlit

# Check status
sudo systemctl status streamlit

# Should show: Active: active (running)
```

#### Step 6.3: Service Management Commands

```bash
# Check if service is running
sudo systemctl status streamlit

# View real-time logs
sudo journalctl -u streamlit -f

# Stop service
sudo systemctl stop streamlit

# Start service
sudo systemctl start streamlit

# Restart service
sudo systemctl restart streamlit

# Disable auto-start on boot
sudo systemctl disable streamlit
```

---

## 🔧 Part 7: Maintenance and Updates

### Update Your Code

```bash
# SSH into EC2
ssh -i ~/path/to/your-key.pem ec2-user@your-ec2-ip

# Navigate to project
cd ~/Enrich_KEs

# Pull latest changes
git pull origin main

# Activate environment
conda activate enrich_kes

# Update dependencies if needed
conda env update -f environment.yml

# Restart service
sudo systemctl restart streamlit

# Check logs
sudo journalctl -u streamlit -f
```

### Update Python Packages

```bash
cd ~/Enrich_KEs
conda activate enrich_kes

# Update specific package
pip install --upgrade package-name

# Or update from requirements.txt
pip install -r requirements.txt --upgrade

# Restart service
sudo systemctl restart streamlit
```

### Add Bioinformatics Tools

```bash
# Edit environment.yml
cd ~/Enrich_KEs
nano environment.yml

# Uncomment or add tools, for example:
#   - fastqc
#   - samtools

# Update environment
conda activate enrich_kes
conda env update -f environment.yml

# Verify installation
fastqc --version

# Restart app
sudo systemctl restart streamlit
```

### Backup Your Data

```bash
# Create backup directory
mkdir -p ~/backups

# Backup data directory
cd ~/Enrich_KEs
tar -czf ~/backups/enrich_kes_data_$(date +%Y%m%d).tar.gz data/

# List backups
ls -lh ~/backups/
```

---

## Part 8: Troubleshooting

### App Won't Start

```bash
# Check service status
sudo systemctl status streamlit

# View detailed logs
sudo journalctl -u streamlit -n 100 --no-pager

# Check if Python can import modules
conda activate enrich_kes
python -c "import streamlit, pandas, plotly"

# Manually run to see errors
cd ~/Enrich_KEs
conda activate enrich_kes
streamlit run app.py --server.port 8501 --server.address 0.0.0.0
```

### Can't Access from Browser

```bash
# 1. Check if app is running locally
curl http://localhost:8501

# 2. Check if port 8501 is open
sudo netstat -tlnp | grep 8501

# 3. Verify EC2 security group (see Part 4)

# 4. Get your public IP
curl http://checkip.amazonaws.com

# 5. Try accessing: http://YOUR-IP:8501
```

### Service Keeps Crashing

```bash
# Check logs for errors
sudo journalctl -u streamlit -n 50

# Common issues:
# - Missing data files
# - Wrong Python path
# - Module import errors
# - Port already in use

# Kill any processes using port 8501
sudo kill $(sudo lsof -t -i:8501)

# Restart service
sudo systemctl restart streamlit
```

### Conda Environment Issues

```bash
# List all environments
conda env list

# Recreate environment
conda deactivate
conda env remove -n enrich_kes
conda env create -f environment.yml

# Update service and restart
sudo systemctl restart streamlit
```

### Memory Issues (EC2 Too Small)

```bash
# Check memory usage
free -h

# Check disk usage
df -h

# If running out of memory:
# - Upgrade to larger EC2 instance (t3.medium or larger)
# - Add swap space:
sudo dd if=/dev/zero of=/swapfile bs=1M count=2048
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

---

## Part 9: Monitoring

### Check App Health

```bash
# Service status
sudo systemctl status streamlit

# CPU and memory usage
top
# Press 'q' to quit

# App logs (live)
sudo journalctl -u streamlit -f

# Check specific errors
sudo journalctl -u streamlit | grep -i error

# App access logs
tail -f ~/.streamlit/logs/
```

### Performance Monitoring

```bash
# Install monitoring tools
sudo yum install htop -y

# Monitor system resources
htop

# Check network connections
ss -tulpn | grep 8501

# Disk I/O
iostat -x 1
```

---

## Part 10: Complete Deployment Script

Save this as `deploy.sh` for automated deployment:

```bash
#!/bin/bash
set -e

echo "============================================"
echo "  Enrich_KEs Deployment Script (Conda)"
echo "============================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

# Check if running as ec2-user
if [ "$USER" != "ec2-user" ]; then
    print_error "Please run this script as ec2-user"
    exit 1
fi

# Update system
print_status "Updating system packages..."
sudo yum update -y
sudo yum install git wget -y

# Install Miniconda if not present
if [ ! -d "$HOME/miniconda3" ]; then
    print_status "Installing Miniconda..."
    cd ~
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    ~/miniconda3/bin/conda init bash
    source ~/.bashrc
    print_status "Miniconda installed successfully"
else
    print_status "Miniconda already installed"
fi

# Add conda to path for this script
export PATH="$HOME/miniconda3/bin:$PATH"

# Configure conda channels
print_status "Configuring conda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Clone or update repository
cd ~
if [ -d "Enrich_KEs" ]; then
    print_status "Updating repository..."
    cd Enrich_KEs
    git pull origin main
else
    print_status "Cloning repository..."
    git clone https://github.com/mei-franzi/Enrich_KEs.git
    cd Enrich_KEs
fi

# Check if environment.yml exists
if [ ! -f "environment.yml" ]; then
    print_warning "environment.yml not found, creating from requirements.txt..."
    cat > environment.yml <<EOF
name: enrich_kes
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.10
  - pip
  - pip:
      - streamlit>=1.28.0
      - pandas>=2.0.0
      - numpy>=1.24.0
      - scipy>=1.10.0
      - statsmodels>=0.14.0
      - openpyxl>=3.1.0
      - plotly>=5.14.0
      - gprofiler-official>=1.0.0
      - matplotlib>=3.7.0
      - pytest>=7.4.0
      - pytest-cov>=4.1.0
EOF
fi

# Create or update conda environment
if conda env list | grep -q "enrich_kes"; then
    print_status "Updating existing conda environment..."
    conda env update -f environment.yml --prune
else
    print_status "Creating new conda environment (this may take 5-10 minutes)..."
    conda env create -f environment.yml
fi

# Verify installation
print_status "Verifying installation..."
source $HOME/miniconda3/bin/activate enrich_kes
python -c "import streamlit, pandas, plotly" && print_status "All Python packages imported successfully"

# Create systemd service
print_status "Creating systemd service..."
sudo tee /etc/systemd/system/streamlit.service > /dev/null <<EOF
[Unit]
Description=Streamlit Enrich_KEs App
After=network.target

[Service]
Type=simple
User=ec2-user
WorkingDirectory=/home/ec2-user/Enrich_KEs
Environment="PATH=/home/ec2-user/miniconda3/envs/enrich_kes/bin:/home/ec2-user/miniconda3/bin:/usr/local/bin:/usr/bin:/bin"
ExecStart=/home/ec2-user/miniconda3/envs/enrich_kes/bin/streamlit run app.py --server.port 8501 --server.address 0.0.0.0
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

# Reload systemd
print_status "Reloading systemd..."
sudo systemctl daemon-reload

# Enable and start service
print_status "Enabling and starting Streamlit service..."
sudo systemctl enable streamlit
sudo systemctl restart streamlit

# Wait a moment for service to start
sleep 3

# Check status
if sudo systemctl is-active --quiet streamlit; then
    print_status "Service is running!"
else
    print_error "Service failed to start. Check logs with: sudo journalctl -u streamlit -n 50"
    exit 1
fi

# Get public IP
PUBLIC_IP=$(curl -s http://checkip.amazonaws.com)

echo ""
echo "============================================"
echo "  Deployment Complete! 🎉"
echo "============================================"
echo ""
echo "📱 Access your app at: http://$PUBLIC_IP:8501"
echo ""
echo "Useful commands:"
echo "  Check status: sudo systemctl status streamlit"
echo "  View logs:    sudo journalctl -u streamlit -f"
echo "  Restart:      sudo systemctl restart streamlit"
echo "  Stop:         sudo systemctl stop streamlit"
echo ""
echo "To update your code:"
echo "  cd ~/Enrich_KEs"
echo "  git pull origin main"
echo "  sudo systemctl restart streamlit"
echo ""
print_status "Don't forget to configure your Security Group for port 8501!"
```

**Make it executable and run:**

```bash
chmod +x deploy.sh
./deploy.sh
```

---

## 📝 Quick Reference

### Daily Commands

```bash
# Check app status
sudo systemctl status streamlit

# View logs
sudo journalctl -u streamlit -f

# Restart app
sudo systemctl restart streamlit

# Update code
cd ~/Enrich_KEs && git pull && sudo systemctl restart streamlit
```

### Conda Commands

```bash
# List environments
conda env list

# Activate environment
conda activate enrich_kes

# List packages
conda list

# Update packages
conda env update -f environment.yml
```

### Git Commands

```bash
# Check for updates
git fetch origin main

# Pull updates
git pull origin main

# Check current branch
git branch

# View recent changes
git log --oneline -5
```

---

## Next Steps

1. ✅ Follow this guide to deploy your app
2. 📱 Test the app thoroughly at `http://your-ip:8501`
3. 🔒 Configure HTTPS (optional, see bonus section below)
4. 🌐 Set up custom domain (optional)
5. 📊 Set up monitoring and backups
6. 👥 Share with your collaborator

---

## Bonus: Advanced Configurations

### Enable HTTPS with Let's Encrypt

```bash
# Install certbot
sudo yum install certbot -y

# Get certificate (requires domain name)
sudo certbot certonly --standalone -d yourdomain.com

# Configure Nginx as reverse proxy (see full tutorial online)
```

### Set Up Custom Domain

1. Buy domain (Route 53, Namecheap, etc.)
2. Point A record to EC2 IP
3. Set up Nginx reverse proxy
4. Enable HTTPS with certbot

### Auto-Backup Data

```bash
# Create backup script
cat > ~/backup.sh <<'EOF'
#!/bin/bash
cd ~/Enrich_KEs
tar -czf ~/backups/data_$(date +%Y%m%d_%H%M%S).tar.gz data/
# Keep only last 7 backups
ls -t ~/backups/*.tar.gz | tail -n +8 | xargs -r rm
EOF

chmod +x ~/backup.sh

# Add to crontab (daily at 2 AM)
(crontab -l 2>/dev/null; echo "0 2 * * * ~/backup.sh") | crontab -
```

---

## Checklist

Before going live, verify:

- [ ] Conda environment created and activated
- [ ] All packages installed successfully
- [ ] App runs locally on EC2
- [ ] Security group allows port 8501
- [ ] App accessible from browser
- [ ] Systemd service running
- [ ] Service auto-starts on reboot
- [ ] Logs are clean (no errors)
- [ ] All features work (upload, enrichment, etc.)
- [ ] Data files are present in `data/` directory

---

## Support

If you encounter issues:

1. Check logs: `sudo journalctl -u streamlit -n 100`
2. Test manually: `conda activate enrich_kes && streamlit run app.py`
3. Review troubleshooting section
4. Check AWS security group settings
5. Verify all files are present in `~/Enrich_KEs`

---

**Last Updated:** 2025-01-07  
**Tested On:** Amazon Linux 2, EC2 t3.small/medium instances

