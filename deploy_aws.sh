#!/bin/bash
set -e

echo "============================================"
echo "  Enrich_KEs AWS Deployment Script"
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
print_warning "Don't forget to configure your Security Group for port 8501!"
echo ""

