# AWS EC2 Deployment Steps - What We Did

## Instance Details
- Instance Type: t2.micro
- OS: Amazon Linux 2023
- Storage: 8 GB EBS volume
- RAM: 1 GB

## Step 1: Connect to EC2

```bash
ssh -i ~/path/to/your-key.pem ec2-user@your-ec2-ip
```

## Step 2: Install Git

```bash
sudo yum install git -y
```

## Step 3: Clone Repository

For private repository, use GitHub Personal Access Token:

```bash
cd ~
git clone https://github.com/mei-franzi/Enrich_KEs.git
# When prompted:
# Username: mei-franzi
# Password: [paste your GitHub Personal Access Token]
cd Enrich_KEs
```

## Step 4: Add Swap Space (for memory)

```bash
sudo dd if=/dev/zero of=/swapfile bs=1M count=4096
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

Verify:
```bash
free -h
```

## Step 5: Attempted Conda Installation (Failed)

```bash
# Downloaded and installed Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
~/miniconda3/bin/conda init bash
source ~/.bashrc

# Tried to create environment
conda env create -f environment.yml
```

**Result:** Failed due to disk space. The 8 GB storage filled to 100%.

## Step 6: Clean Up Disk Space

```bash
# Remove failed conda installation
rm -rf ~/miniconda3

# Clear system caches
sudo yum clean all
rm -rf ~/.cache

# Verify space available
df -h
# Result: 2.1 GB available (74% used)
```

## Step 7: Use Python venv Instead

```bash
cd ~/Enrich_KEs

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

Installation completed successfully in approximately 2-3 minutes.

## Step 8: Run the Application

```bash
streamlit run app.py --server.port 8501 --server.address 0.0.0.0
```

## Step 9: Configure Security Group

In AWS Console:
1. Go to EC2 Dashboard
2. Select your instance
3. Security tab > Security Groups > Edit inbound rules
4. Add rule:
   - Type: Custom TCP
   - Port: 8501
   - Source: 0.0.0.0/0 (or specific IP addresses)

## Access the Application

```
http://your-ec2-public-ip:8501
```

## Key Issues and Solutions

### Issue 1: Out of Memory during Conda Installation
**Solution:** Added 4 GB swap space

### Issue 2: Out of Disk Space (100% full)
**Cause:** Conda installation requires ~3 GB, but only 8 GB total disk
**Solution:** Removed conda, used pip/venv instead (~500 MB)

### Issue 3: Corrupted venv
**Cause:** Disk full during venv creation
**Solution:** Deleted and recreated venv after freeing disk space

## Final Working Configuration

- Python: 3.13 (system Python)
- Package Manager: pip
- Virtual Environment: venv
- Dependencies: Installed from requirements.txt
- Disk Usage: ~6 GB used, 2 GB free
- Memory Usage: 1 GB RAM + 4 GB swap

## Recommendations for Future

**For Production or Team Use:**
- Upgrade to t3.small (2 GB RAM, 20+ GB storage)
- This would allow conda installation for bioinformatics tools
- Current venv setup works fine for current application needs

